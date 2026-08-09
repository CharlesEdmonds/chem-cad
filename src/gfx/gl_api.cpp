#include "gfx/gl_api.hpp"

#include <string>
#include <type_traits>

namespace chemcad::gfx {

#define CHEMCAD_GL_PROCS(X)                       \
  X(glGetString)                                  \
  X(glGetError)                                   \
  X(glGetIntegerv)                                \
  X(glGetBooleanv)                                \
  X(glGetFloatv)                                  \
  X(glIsEnabled)                                  \
  X(glEnable)                                     \
  X(glDisable)                                    \
  X(glViewport)                                   \
  X(glClear)                                      \
  X(glClearBufferfv)                              \
  X(glClearColor)                                 \
  X(glColorMask)                                  \
  X(glDepthFunc)                                  \
  X(glDepthMask)                                  \
  X(glBlendFuncSeparate)                          \
  X(glBlendEquationSeparate)                      \
  X(glCullFace)                                   \
  X(glFrontFace)                                  \
  X(glGenBuffers)                                 \
  X(glDeleteBuffers)                              \
  X(glBindBuffer)                                 \
  X(glBufferData)                                 \
  X(glGenVertexArrays)                            \
  X(glDeleteVertexArrays)                         \
  X(glBindVertexArray)                            \
  X(glEnableVertexAttribArray)                    \
  X(glVertexAttribPointer)                        \
  X(glVertexAttribDivisor)                        \
  X(glCreateShader)                               \
  X(glShaderSource)                               \
  X(glCompileShader)                              \
  X(glGetShaderiv)                                \
  X(glGetShaderInfoLog)                           \
  X(glDeleteShader)                               \
  X(glCreateProgram)                              \
  X(glAttachShader)                               \
  X(glLinkProgram)                                \
  X(glGetProgramiv)                               \
  X(glGetProgramInfoLog)                          \
  X(glDeleteProgram)                              \
  X(glUseProgram)                                 \
  X(glGetUniformLocation)                         \
  X(glUniform1i)                                  \
  X(glUniform1f)                                  \
  X(glUniform2f)                                  \
  X(glUniform3f)                                  \
  X(glUniform4f)                                  \
  X(glUniform1fv)                                 \
  X(glUniformMatrix3fv)                           \
  X(glUniformMatrix4fv)                           \
  X(glGenTextures)                                \
  X(glDeleteTextures)                             \
  X(glActiveTexture)                              \
  X(glBindTexture)                                \
  X(glTexImage2D)                                 \
  X(glTexParameteri)                              \
  X(glGenFramebuffers)                            \
  X(glDeleteFramebuffers)                         \
  X(glBindFramebuffer)                            \
  X(glFramebufferTexture2D)                       \
  X(glCheckFramebufferStatus)                     \
  X(glDrawBuffers)                                \
  X(glGenRenderbuffers)                           \
  X(glDeleteRenderbuffers)                        \
  X(glBindRenderbuffer)                           \
  X(glRenderbufferStorage)                        \
  X(glFramebufferRenderbuffer)                    \
  X(glDrawArrays)                                 \
  X(glDrawArraysInstanced)                        \
  X(glDrawElements)

#define CHEMCAD_DEFINE_PROC(name) decltype(name) name = nullptr;
CHEMCAD_GL_PROCS(CHEMCAD_DEFINE_PROC)
#undef CHEMCAD_DEFINE_PROC

namespace {

bool loaded = false;
std::string versionString;
std::string rendererString;

void clearProcs() {
#define CHEMCAD_CLEAR_PROC(name) name = nullptr;
  CHEMCAD_GL_PROCS(CHEMCAD_CLEAR_PROC)
#undef CHEMCAD_CLEAR_PROC
  loaded = false;
  versionString.clear();
  rendererString.clear();
}

struct ProcTable {
#define CHEMCAD_TABLE_PROC(name) decltype(::chemcad::gfx::name) name = nullptr;
  CHEMCAD_GL_PROCS(CHEMCAD_TABLE_PROC)
#undef CHEMCAD_TABLE_PROC
};

}  // namespace

bool loadGl(GlProcLoader loader) {
  clearProcs();
  if (loader == nullptr) return false;

  ProcTable table;
  bool complete = true;
  auto resolve = [&](auto& slot, const char* name) {
    using Proc = std::remove_reference_t<decltype(slot)>;
    slot = reinterpret_cast<Proc>(loader(name));
    if (slot == nullptr) complete = false;
  };

#define CHEMCAD_LOAD_PROC(name) resolve(table.name, #name);
  CHEMCAD_GL_PROCS(CHEMCAD_LOAD_PROC)
#undef CHEMCAD_LOAD_PROC

  if (!complete) return false;

#define CHEMCAD_ASSIGN_PROC(name) name = table.name;
  CHEMCAD_GL_PROCS(CHEMCAD_ASSIGN_PROC)
#undef CHEMCAD_ASSIGN_PROC

  const GLubyte* version = glGetString(GL_VERSION);
  const GLubyte* renderer = glGetString(GL_RENDERER);
  versionString = version != nullptr ? reinterpret_cast<const char*>(version) : "Unavailable";
  rendererString = renderer != nullptr ? reinterpret_cast<const char*>(renderer) : "Unavailable";
  loaded = true;
  return true;
}

bool glLoaded() { return loaded; }

const char* glVersionString() { return versionString.empty() ? "Unavailable" : versionString.c_str(); }

const char* glRendererString() {
  return rendererString.empty() ? "Unavailable" : rendererString.c_str();
}

#undef CHEMCAD_GL_PROCS

}  // namespace chemcad::gfx
