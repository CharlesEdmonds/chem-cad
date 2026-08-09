#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

#include "gfx/gl_api.hpp"

namespace chemcad::gfx {

struct ComputeCapabilities {
  bool available = false;
  int majorVersion = 0;
  int minorVersion = 0;
  int maxWorkGroupInvocations = 0;
  int maxStorageBlocks = 0;
  int maxStorageBufferBindings = 0;
  std::int64_t maxStorageBlockBytes = 0;
  std::string reason;
};

// Safe before GL initialisation and without a current context. A false result
// means callers must retain their non-compute path.
ComputeCapabilities probeCompute();

class ComputeProgram {
 public:
  ComputeProgram() = default;
  ~ComputeProgram();
  ComputeProgram(const ComputeProgram&) = delete;
  ComputeProgram& operator=(const ComputeProgram&) = delete;
  ComputeProgram(ComputeProgram&& other) noexcept;
  ComputeProgram& operator=(ComputeProgram&& other) noexcept;

  bool compile(std::string_view source, std::string* errorLog = nullptr);
  void reset();
  bool valid() const { return program_ != 0; }
  GLuint id() const { return program_; }

  void use() const;
  void dispatch(GLuint x, GLuint y = 1, GLuint z = 1) const;
  void setInt(const char* name, int value) const;
  void setUInt(const char* name, std::uint32_t value) const;
  void setFloat(const char* name, float value) const;
  void setVec3(const char* name, float x, float y, float z) const;
  void setIVec3(const char* name, int x, int y, int z) const;
  void bindStorageBlock(GLuint blockIndex, GLuint binding) const;

  static void memoryBarrier(
      GLbitfield bits = GL_SHADER_STORAGE_BARRIER_BIT | GL_BUFFER_UPDATE_BARRIER_BIT);

 private:
  GLint uniform(const char* name) const;
  GLuint program_ = 0;
};

class StorageBuffer {
 public:
  StorageBuffer() = default;
  ~StorageBuffer();
  StorageBuffer(const StorageBuffer&) = delete;
  StorageBuffer& operator=(const StorageBuffer&) = delete;
  StorageBuffer(StorageBuffer&& other) noexcept;
  StorageBuffer& operator=(StorageBuffer&& other) noexcept;

  bool resize(std::size_t bytes, GLenum usage = GL_DYNAMIC_DRAW);
  bool upload(const void* data, std::size_t bytes, std::size_t offset = 0);
  bool download(void* data, std::size_t bytes, std::size_t offset = 0) const;
  bool clearUInt(std::uint32_t value = 0);
  void bind(GLuint binding) const;
  void reset();

  bool valid() const { return buffer_ != 0; }
  GLuint id() const { return buffer_; }
  std::size_t size() const { return bytes_; }

 private:
  GLuint buffer_ = 0;
  std::size_t bytes_ = 0;
};

}  // namespace chemcad::gfx
