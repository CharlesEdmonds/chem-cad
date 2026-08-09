#include "gfx/compute.hpp"

#include <algorithm>
#include <cstring>
#include <utility>
#include <vector>

namespace chemcad::gfx {
namespace {

constexpr int kRequiredStorageBindings = 12;
constexpr int kRequiredInvocations = 128;
constexpr std::int64_t kRequiredStorageBlockBytes = 16 * 1024 * 1024;

bool computeEntryPointsPresent() {
  return glGetIntegeri_v != nullptr && glGetInteger64v != nullptr &&
         glBufferSubData != nullptr && glGetBufferSubData != nullptr &&
         glBindBufferBase != nullptr && glMapBufferRange != nullptr &&
         glUnmapBuffer != nullptr && glClearBufferData != nullptr &&
         glUniform1ui != nullptr && glUniform3i != nullptr &&
         glDispatchCompute != nullptr && glMemoryBarrier != nullptr &&
         glShaderStorageBlockBinding != nullptr;
}

std::string shaderLog(GLuint shader) {
  GLint length = 0;
  glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
  if (length <= 1) return {};
  std::vector<GLchar> bytes(static_cast<std::size_t>(length));
  GLsizei written = 0;
  glGetShaderInfoLog(shader, length, &written, bytes.data());
  return std::string(bytes.data(), static_cast<std::size_t>(std::max(0, written)));
}

std::string programLog(GLuint program) {
  GLint length = 0;
  glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);
  if (length <= 1) return {};
  std::vector<GLchar> bytes(static_cast<std::size_t>(length));
  GLsizei written = 0;
  glGetProgramInfoLog(program, length, &written, bytes.data());
  return std::string(bytes.data(), static_cast<std::size_t>(std::max(0, written)));
}

}  // namespace

ComputeCapabilities probeCompute() {
  ComputeCapabilities result;
  if (!glLoaded() || glGetString == nullptr || glGetIntegerv == nullptr) {
    result.reason = "OpenGL is not loaded or no context is current";
    return result;
  }
  if (!computeEntryPointsPresent()) {
    result.reason = "required OpenGL compute entry points are unavailable";
    return result;
  }
  if (glGetString(GL_VERSION) == nullptr) {
    result.reason = "no current OpenGL context";
    return result;
  }

  glGetIntegerv(GL_MAJOR_VERSION, &result.majorVersion);
  glGetIntegerv(GL_MINOR_VERSION, &result.minorVersion);
  if (result.majorVersion < 4 ||
      (result.majorVersion == 4 && result.minorVersion < 3)) {
    result.reason = "OpenGL 4.3 or newer is required";
    return result;
  }

  glGetIntegerv(GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS,
                &result.maxWorkGroupInvocations);
  glGetIntegerv(GL_MAX_COMPUTE_SHADER_STORAGE_BLOCKS,
                &result.maxStorageBlocks);
  glGetIntegerv(GL_MAX_SHADER_STORAGE_BUFFER_BINDINGS,
                &result.maxStorageBufferBindings);
  glGetInteger64v(GL_MAX_SHADER_STORAGE_BLOCK_SIZE,
                  &result.maxStorageBlockBytes);

  GLint groupCount = 0;
  GLint groupSize = 0;
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &groupCount);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 0, &groupSize);
  if (result.maxWorkGroupInvocations < kRequiredInvocations || groupSize < 128 ||
      groupCount < 1) {
    result.reason = "compute work-group limits are below the solver minimum";
    return result;
  }
  if (result.maxStorageBlocks < kRequiredStorageBindings ||
      result.maxStorageBufferBindings < kRequiredStorageBindings) {
    result.reason = "too few shader-storage buffer bindings";
    return result;
  }
  if (result.maxStorageBlockBytes < kRequiredStorageBlockBytes) {
    result.reason = "shader-storage blocks are smaller than 16 MiB";
    return result;
  }

  result.available = true;
  return result;
}

ComputeProgram::~ComputeProgram() { reset(); }

ComputeProgram::ComputeProgram(ComputeProgram&& other) noexcept
    : program_(std::exchange(other.program_, 0)) {}

ComputeProgram& ComputeProgram::operator=(ComputeProgram&& other) noexcept {
  if (this != &other) {
    reset();
    program_ = std::exchange(other.program_, 0);
  }
  return *this;
}

bool ComputeProgram::compile(std::string_view source, std::string* errorLog) {
  reset();
  if (!probeCompute().available || source.empty()) {
    if (errorLog != nullptr) *errorLog = "OpenGL compute is unavailable";
    return false;
  }

  const GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
  if (shader == 0) {
    if (errorLog != nullptr) *errorLog = "glCreateShader returned zero";
    return false;
  }
  const GLchar* text = source.data();
  const GLint length = static_cast<GLint>(source.size());
  glShaderSource(shader, 1, &text, &length);
  glCompileShader(shader);
  GLint compiled = GL_FALSE;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
  if (compiled != GL_TRUE) {
    if (errorLog != nullptr) *errorLog = shaderLog(shader);
    glDeleteShader(shader);
    return false;
  }

  const GLuint program = glCreateProgram();
  if (program == 0) {
    if (errorLog != nullptr) *errorLog = "glCreateProgram returned zero";
    glDeleteShader(shader);
    return false;
  }
  glAttachShader(program, shader);
  glLinkProgram(program);
  glDeleteShader(shader);
  GLint linked = GL_FALSE;
  glGetProgramiv(program, GL_LINK_STATUS, &linked);
  if (linked != GL_TRUE) {
    if (errorLog != nullptr) *errorLog = programLog(program);
    glDeleteProgram(program);
    return false;
  }

  program_ = program;
  if (errorLog != nullptr) errorLog->clear();
  return true;
}

void ComputeProgram::reset() {
  if (program_ != 0 && glDeleteProgram != nullptr) glDeleteProgram(program_);
  program_ = 0;
}

void ComputeProgram::use() const {
  if (program_ != 0 && glUseProgram != nullptr) glUseProgram(program_);
}

void ComputeProgram::dispatch(GLuint x, GLuint y, GLuint z) const {
  if (program_ == 0 || glDispatchCompute == nullptr || x == 0 || y == 0 || z == 0) return;
  use();
  glDispatchCompute(x, y, z);
}

GLint ComputeProgram::uniform(const char* name) const {
  return program_ != 0 && name != nullptr && glGetUniformLocation != nullptr
             ? glGetUniformLocation(program_, name)
             : -1;
}

void ComputeProgram::setInt(const char* name, int value) const {
  const GLint location = uniform(name);
  if (location >= 0) glUniform1i(location, value);
}

void ComputeProgram::setUInt(const char* name, std::uint32_t value) const {
  const GLint location = uniform(name);
  if (location >= 0) glUniform1ui(location, value);
}

void ComputeProgram::setFloat(const char* name, float value) const {
  const GLint location = uniform(name);
  if (location >= 0) glUniform1f(location, value);
}

void ComputeProgram::setVec3(const char* name, float x, float y, float z) const {
  const GLint location = uniform(name);
  if (location >= 0) glUniform3f(location, x, y, z);
}

void ComputeProgram::setIVec3(const char* name, int x, int y, int z) const {
  const GLint location = uniform(name);
  if (location >= 0) glUniform3i(location, x, y, z);
}

void ComputeProgram::bindStorageBlock(GLuint blockIndex, GLuint binding) const {
  if (program_ != 0 && glShaderStorageBlockBinding != nullptr) {
    glShaderStorageBlockBinding(program_, blockIndex, binding);
  }
}

void ComputeProgram::memoryBarrier(GLbitfield bits) {
  if (glMemoryBarrier != nullptr) glMemoryBarrier(bits);
}

StorageBuffer::~StorageBuffer() { reset(); }

StorageBuffer::StorageBuffer(StorageBuffer&& other) noexcept
    : buffer_(std::exchange(other.buffer_, 0)),
      bytes_(std::exchange(other.bytes_, 0)) {}

StorageBuffer& StorageBuffer::operator=(StorageBuffer&& other) noexcept {
  if (this != &other) {
    reset();
    buffer_ = std::exchange(other.buffer_, 0);
    bytes_ = std::exchange(other.bytes_, 0);
  }
  return *this;
}

bool StorageBuffer::resize(std::size_t bytes, GLenum usage) {
  const ComputeCapabilities capabilities = probeCompute();
  if (!capabilities.available || bytes == 0 ||
      bytes > static_cast<std::size_t>(capabilities.maxStorageBlockBytes)) {
    return false;
  }
  if (buffer_ != 0 && bytes_ == bytes) return true;
  if (buffer_ == 0) glGenBuffers(1, &buffer_);
  if (buffer_ == 0) return false;
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffer_);
  glBufferData(GL_SHADER_STORAGE_BUFFER, static_cast<GLsizeiptr>(bytes), nullptr, usage);
  bytes_ = bytes;
  return glGetError == nullptr || glGetError() == GL_NO_ERROR;
}

bool StorageBuffer::upload(const void* data, std::size_t bytes, std::size_t offset) {
  if (buffer_ == 0 || data == nullptr || offset > bytes_ || bytes > bytes_ - offset) return false;
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffer_);
  glBufferSubData(GL_SHADER_STORAGE_BUFFER, static_cast<GLintptr>(offset),
                  static_cast<GLsizeiptr>(bytes), data);
  return glGetError == nullptr || glGetError() == GL_NO_ERROR;
}

bool StorageBuffer::download(void* data, std::size_t bytes, std::size_t offset) const {
  if (buffer_ == 0 || data == nullptr || offset > bytes_ || bytes > bytes_ - offset) return false;
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffer_);
  void* mapped = glMapBufferRange(GL_SHADER_STORAGE_BUFFER,
                                  static_cast<GLintptr>(offset),
                                  static_cast<GLsizeiptr>(bytes), GL_MAP_READ_BIT);
  if (mapped == nullptr) return false;
  std::memcpy(data, mapped, bytes);
  return glUnmapBuffer(GL_SHADER_STORAGE_BUFFER) == GL_TRUE;
}

bool StorageBuffer::clearUInt(std::uint32_t value) {
  if (buffer_ == 0) return false;
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffer_);
  glClearBufferData(GL_SHADER_STORAGE_BUFFER, GL_R32UI, GL_RED_INTEGER,
                    GL_UNSIGNED_INT, &value);
  return glGetError == nullptr || glGetError() == GL_NO_ERROR;
}

void StorageBuffer::bind(GLuint binding) const {
  if (buffer_ != 0 && glBindBufferBase != nullptr) {
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, binding, buffer_);
  }
}

void StorageBuffer::reset() {
  if (buffer_ != 0 && glDeleteBuffers != nullptr) glDeleteBuffers(1, &buffer_);
  buffer_ = 0;
  bytes_ = 0;
}

}  // namespace chemcad::gfx
