#pragma once
// Minimal OpenGL 3.3 core function loader.
//
// The application creates a GL 3.3 core context (src/app/main.cpp) and Dear
// ImGui's backend carries its own private loader, whose header states plainly
// that the rest of the application must use a different one. Rather than adding
// a code generator or a whole loader library for the ~30 entry points the fluid
// renderer needs, this declares exactly those and resolves them through the
// context's own getProcAddress.
//
// Everything here is a no-op until `loadGl` succeeds, and every pointer is null
// until then, so a headless build that never calls it cannot accidentally issue
// a GL call.

#include <cstddef>
#include <cstdint>

namespace chemcad::gfx {
using GLenum = std::uint32_t;
using GLboolean = std::uint8_t;
using GLbitfield = std::uint32_t;
using GLint = int;
using GLsizei = int;
using GLuint = std::uint32_t;
using GLfloat = float;
using GLchar = char;
using GLubyte = unsigned char;
using GLintptr = std::intptr_t;
using GLsizeiptr = std::intptr_t;


// The system's GL 1.1 header (pulled in by GLFW, or directly by any translation
// unit that needs glReadPixels) defines the classic enums as MACROS, which
// ignore namespaces and would break the typed constants below. Undefining them
// here makes this header safe to include anywhere, in any order; the typed
// constants that replace them live in chemcad::gfx and cannot collide again.
#undef GL_ACTIVE_TEXTURE
#undef GL_ARRAY_BUFFER
#undef GL_ARRAY_BUFFER_BINDING
#undef GL_BACK
#undef GL_BLEND
#undef GL_BLEND_DST_ALPHA
#undef GL_BLEND_DST_RGB
#undef GL_BLEND_EQUATION_ALPHA
#undef GL_BLEND_EQUATION_RGB
#undef GL_BLEND_SRC_ALPHA
#undef GL_BLEND_SRC_RGB
#undef GL_CCW
#undef GL_CLAMP_TO_EDGE
#undef GL_COLOR
#undef GL_COLOR_ATTACHMENT0
#undef GL_COLOR_ATTACHMENT1
#undef GL_COLOR_BUFFER_BIT
#undef GL_COLOR_CLEAR_VALUE
#undef GL_COLOR_WRITEMASK
#undef GL_COMPILE_STATUS
#undef GL_CULL_FACE
#undef GL_CULL_FACE_MODE
#undef GL_CURRENT_PROGRAM
#undef GL_DEPTH
#undef GL_DEPTH24_STENCIL8
#undef GL_DEPTH_BUFFER_BIT
#undef GL_DEPTH_FUNC
#undef GL_DEPTH_STENCIL_ATTACHMENT
#undef GL_DEPTH_TEST
#undef GL_DEPTH_WRITEMASK
#undef GL_DRAW_FRAMEBUFFER
#undef GL_DRAW_FRAMEBUFFER_BINDING
#undef GL_ELEMENT_ARRAY_BUFFER
#undef GL_FALSE
#undef GL_FLOAT
#undef GL_FRAGMENT_SHADER
#undef GL_FRAMEBUFFER
#undef GL_FRAMEBUFFER_COMPLETE
#undef GL_FRAMEBUFFER_SRGB
#undef GL_FRONT
#undef GL_FRONT_FACE
#undef GL_FUNC_ADD
#undef GL_INFO_LOG_LENGTH
#undef GL_LESS
#undef GL_LINEAR
#undef GL_LINK_STATUS
#undef GL_NEAREST
#undef GL_NO_ERROR
#undef GL_ONE
#undef GL_ONE_MINUS_SRC_ALPHA
#undef GL_R16F
#undef GL_R32F
#undef GL_READ_FRAMEBUFFER
#undef GL_READ_FRAMEBUFFER_BINDING
#undef GL_RED
#undef GL_RENDERBUFFER
#undef GL_RENDERBUFFER_BINDING
#undef GL_RENDERER
#undef GL_RG
#undef GL_RG16F
#undef GL_RGBA
#undef GL_RGBA8
#undef GL_SCISSOR_TEST
#undef GL_SRC_ALPHA
#undef GL_STATIC_DRAW
#undef GL_STENCIL_BUFFER_BIT
#undef GL_STENCIL_TEST
#undef GL_STREAM_DRAW
#undef GL_TEXTURE0
#undef GL_TEXTURE_2D
#undef GL_TEXTURE_BINDING_2D
#undef GL_TEXTURE_MAG_FILTER
#undef GL_TEXTURE_MIN_FILTER
#undef GL_TEXTURE_WRAP_S
#undef GL_TEXTURE_WRAP_T
#undef GL_TRIANGLES
#undef GL_TRUE
#undef GL_UNSIGNED_BYTE
#undef GL_UNSIGNED_INT
#undef GL_VERSION
#undef GL_VERTEX_ARRAY_BINDING
#undef GL_VERTEX_SHADER
#undef GL_VIEWPORT
#undef GL_ZERO

inline constexpr GLboolean GL_FALSE = 0;
inline constexpr GLboolean GL_TRUE = 1;
inline constexpr GLenum GL_NO_ERROR = 0;
inline constexpr GLenum GL_VERSION = 0x1F02;
inline constexpr GLenum GL_RENDERER = 0x1F01;
inline constexpr GLenum GL_ARRAY_BUFFER = 0x8892;
inline constexpr GLenum GL_ARRAY_BUFFER_BINDING = 0x8894;
inline constexpr GLenum GL_ELEMENT_ARRAY_BUFFER = 0x8893;
inline constexpr GLenum GL_STATIC_DRAW = 0x88E4;
inline constexpr GLenum GL_STREAM_DRAW = 0x88E0;
inline constexpr GLenum GL_FLOAT = 0x1406;
inline constexpr GLenum GL_UNSIGNED_INT = 0x1405;
inline constexpr GLenum GL_TRIANGLES = 0x0004;
inline constexpr GLenum GL_VERTEX_SHADER = 0x8B31;
inline constexpr GLenum GL_FRAGMENT_SHADER = 0x8B30;
inline constexpr GLenum GL_COMPILE_STATUS = 0x8B81;
inline constexpr GLenum GL_LINK_STATUS = 0x8B82;
inline constexpr GLenum GL_INFO_LOG_LENGTH = 0x8B84;
inline constexpr GLenum GL_CURRENT_PROGRAM = 0x8B8D;
inline constexpr GLenum GL_VERTEX_ARRAY_BINDING = 0x85B5;
inline constexpr GLenum GL_TEXTURE0 = 0x84C0;
inline constexpr GLenum GL_ACTIVE_TEXTURE = 0x84E0;
inline constexpr GLenum GL_TEXTURE_2D = 0x0DE1;
inline constexpr GLenum GL_TEXTURE_BINDING_2D = 0x8069;
inline constexpr GLenum GL_TEXTURE_MIN_FILTER = 0x2801;
inline constexpr GLenum GL_TEXTURE_MAG_FILTER = 0x2800;
inline constexpr GLenum GL_TEXTURE_WRAP_S = 0x2802;
inline constexpr GLenum GL_TEXTURE_WRAP_T = 0x2803;
inline constexpr GLenum GL_NEAREST = 0x2600;
inline constexpr GLenum GL_LINEAR = 0x2601;
inline constexpr GLenum GL_CLAMP_TO_EDGE = 0x812F;
inline constexpr GLenum GL_RED = 0x1903;
inline constexpr GLenum GL_RG = 0x8227;
inline constexpr GLenum GL_RGBA = 0x1908;
inline constexpr GLenum GL_R16F = 0x822D;
inline constexpr GLenum GL_R32F = 0x822E;
inline constexpr GLenum GL_RG16F = 0x822F;
inline constexpr GLenum GL_RGBA8 = 0x8058;
inline constexpr GLenum GL_UNSIGNED_BYTE = 0x1401;
inline constexpr GLenum GL_FRAMEBUFFER = 0x8D40;
inline constexpr GLenum GL_READ_FRAMEBUFFER = 0x8CA8;
inline constexpr GLenum GL_DRAW_FRAMEBUFFER = 0x8CA9;
inline constexpr GLenum GL_READ_FRAMEBUFFER_BINDING = 0x8CAA;
inline constexpr GLenum GL_DRAW_FRAMEBUFFER_BINDING = 0x8CA6;
inline constexpr GLenum GL_FRAMEBUFFER_COMPLETE = 0x8CD5;
inline constexpr GLenum GL_COLOR_ATTACHMENT0 = 0x8CE0;
inline constexpr GLenum GL_COLOR_ATTACHMENT1 = 0x8CE1;
inline constexpr GLenum GL_DEPTH_STENCIL_ATTACHMENT = 0x821A;
inline constexpr GLenum GL_RENDERBUFFER = 0x8D41;
inline constexpr GLenum GL_DEPTH24_STENCIL8 = 0x88F0;
inline constexpr GLenum GL_COLOR = 0x1800;
inline constexpr GLenum GL_RENDERBUFFER_BINDING = 0x8CA7;
inline constexpr GLenum GL_DEPTH = 0x1801;
inline constexpr GLbitfield GL_DEPTH_BUFFER_BIT = 0x00000100;
inline constexpr GLbitfield GL_STENCIL_BUFFER_BIT = 0x00000400;
inline constexpr GLenum GL_VIEWPORT = 0x0BA2;
inline constexpr GLbitfield GL_COLOR_BUFFER_BIT = 0x00004000;
inline constexpr GLenum GL_COLOR_CLEAR_VALUE = 0x0C22;
inline constexpr GLenum GL_BLEND = 0x0BE2;
inline constexpr GLenum GL_DEPTH_TEST = 0x0B71;
inline constexpr GLenum GL_SCISSOR_TEST = 0x0C11;
inline constexpr GLenum GL_STENCIL_TEST = 0x0B90;
inline constexpr GLenum GL_FRAMEBUFFER_SRGB = 0x8DB9;
inline constexpr GLenum GL_COLOR_WRITEMASK = 0x0C23;
inline constexpr GLenum GL_CULL_FACE = 0x0B44;
inline constexpr GLenum GL_BLEND_SRC_RGB = 0x80C9;
inline constexpr GLenum GL_BLEND_DST_RGB = 0x80C8;
inline constexpr GLenum GL_BLEND_SRC_ALPHA = 0x80CB;
inline constexpr GLenum GL_BLEND_DST_ALPHA = 0x80CA;
inline constexpr GLenum GL_BLEND_EQUATION_RGB = 0x8009;
inline constexpr GLenum GL_BLEND_EQUATION_ALPHA = 0x883D;
inline constexpr GLenum GL_FUNC_ADD = 0x8006;
inline constexpr GLenum GL_ZERO = 0;
inline constexpr GLenum GL_ONE = 1;
inline constexpr GLenum GL_SRC_ALPHA = 0x0302;
inline constexpr GLenum GL_ONE_MINUS_SRC_ALPHA = 0x0303;
inline constexpr GLenum GL_DEPTH_WRITEMASK = 0x0B72;
inline constexpr GLenum GL_DEPTH_FUNC = 0x0B74;
inline constexpr GLenum GL_LESS = 0x0201;
inline constexpr GLenum GL_CULL_FACE_MODE = 0x0B45;
inline constexpr GLenum GL_FRONT_FACE = 0x0B46;
inline constexpr GLenum GL_FRONT = 0x0404;
inline constexpr GLenum GL_BACK = 0x0405;
inline constexpr GLenum GL_CCW = 0x0901;

#if defined(_WIN32)
#define CHEMCAD_GL_APIENTRY __stdcall
#else
#define CHEMCAD_GL_APIENTRY
#endif

using PFNGLGETSTRINGPROC = const GLubyte* (CHEMCAD_GL_APIENTRY*)(GLenum);
using PFNGLGETERRORPROC = GLenum (CHEMCAD_GL_APIENTRY*)();
using PFNGLGETINTEGERVPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLint*);
using PFNGLGETBOOLEANVPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLboolean*);
using PFNGLGETFLOATVPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLfloat*);
using PFNGLISENABLEDPROC = GLboolean (CHEMCAD_GL_APIENTRY*)(GLenum);
using PFNGLENABLEPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum);
using PFNGLDISABLEPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum);
using PFNGLVIEWPORTPROC = void (CHEMCAD_GL_APIENTRY*)(GLint, GLint, GLsizei, GLsizei);
using PFNGLCLEARPROC = void (CHEMCAD_GL_APIENTRY*)(GLbitfield);
using PFNGLCLEARBUFFERFVPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLint, const GLfloat*);
using PFNGLDEPTHFUNCPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum);
using PFNGLDEPTHMASKPROC = void (CHEMCAD_GL_APIENTRY*)(GLboolean);
using PFNGLBLENDFUNCSEPARATEPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLenum, GLenum, GLenum, GLenum);
using PFNGLBLENDEQUATIONSEPARATEPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLenum);
using PFNGLCULLFACEPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum);
using PFNGLFRONTFACEPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum);
using PFNGLGENBUFFERSPROC = void (CHEMCAD_GL_APIENTRY*)(GLsizei, GLuint*);
using PFNGLDELETEBUFFERSPROC = void (CHEMCAD_GL_APIENTRY*)(GLsizei, const GLuint*);
using PFNGLCOLORMASKPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLboolean, GLboolean, GLboolean, GLboolean);
using PFNGLBINDBUFFERPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLuint);
using PFNGLCLEARCOLORPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLfloat, GLfloat, GLfloat, GLfloat);
using PFNGLBUFFERDATAPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLsizeiptr, const void*, GLenum);
using PFNGLGENVERTEXARRAYSPROC = void (CHEMCAD_GL_APIENTRY*)(GLsizei, GLuint*);
using PFNGLDELETEVERTEXARRAYSPROC = void (CHEMCAD_GL_APIENTRY*)(GLsizei, const GLuint*);
using PFNGLBINDVERTEXARRAYPROC = void (CHEMCAD_GL_APIENTRY*)(GLuint);
using PFNGLENABLEVERTEXATTRIBARRAYPROC = void (CHEMCAD_GL_APIENTRY*)(GLuint);
using PFNGLVERTEXATTRIBPOINTERPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLuint, GLint, GLenum, GLboolean, GLsizei, const void*);
using PFNGLVERTEXATTRIBDIVISORPROC = void (CHEMCAD_GL_APIENTRY*)(GLuint, GLuint);
using PFNGLCREATESHADERPROC = GLuint (CHEMCAD_GL_APIENTRY*)(GLenum);
using PFNGLSHADERSOURCEPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLuint, GLsizei, const GLchar* const*, const GLint*);
using PFNGLCOMPILESHADERPROC = void (CHEMCAD_GL_APIENTRY*)(GLuint);
using PFNGLGETSHADERIVPROC = void (CHEMCAD_GL_APIENTRY*)(GLuint, GLenum, GLint*);
using PFNGLGETSHADERINFOLOGPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLuint, GLsizei, GLsizei*, GLchar*);
using PFNGLDELETESHADERPROC = void (CHEMCAD_GL_APIENTRY*)(GLuint);
using PFNGLCREATEPROGRAMPROC = GLuint (CHEMCAD_GL_APIENTRY*)();
using PFNGLATTACHSHADERPROC = void (CHEMCAD_GL_APIENTRY*)(GLuint, GLuint);
using PFNGLLINKPROGRAMPROC = void (CHEMCAD_GL_APIENTRY*)(GLuint);
using PFNGLGETPROGRAMIVPROC = void (CHEMCAD_GL_APIENTRY*)(GLuint, GLenum, GLint*);
using PFNGLGETPROGRAMINFOLOGPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLuint, GLsizei, GLsizei*, GLchar*);
using PFNGLDELETEPROGRAMPROC = void (CHEMCAD_GL_APIENTRY*)(GLuint);
using PFNGLUSEPROGRAMPROC = void (CHEMCAD_GL_APIENTRY*)(GLuint);
using PFNGLGETUNIFORMLOCATIONPROC = GLint (CHEMCAD_GL_APIENTRY*)(GLuint, const GLchar*);
using PFNGLUNIFORM1IPROC = void (CHEMCAD_GL_APIENTRY*)(GLint, GLint);
using PFNGLUNIFORM1FPROC = void (CHEMCAD_GL_APIENTRY*)(GLint, GLfloat);
using PFNGLUNIFORM2FPROC = void (CHEMCAD_GL_APIENTRY*)(GLint, GLfloat, GLfloat);
using PFNGLUNIFORM3FPROC = void (CHEMCAD_GL_APIENTRY*)(GLint, GLfloat, GLfloat, GLfloat);
using PFNGLUNIFORM4FPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLint, GLfloat, GLfloat, GLfloat, GLfloat);
using PFNGLUNIFORM1FVPROC = void (CHEMCAD_GL_APIENTRY*)(GLint, GLsizei, const GLfloat*);
using PFNGLUNIFORMMATRIX3FVPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLint, GLsizei, GLboolean, const GLfloat*);
using PFNGLUNIFORMMATRIX4FVPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLint, GLsizei, GLboolean, const GLfloat*);
using PFNGLGENTEXTURESPROC = void (CHEMCAD_GL_APIENTRY*)(GLsizei, GLuint*);
using PFNGLDELETETEXTURESPROC = void (CHEMCAD_GL_APIENTRY*)(GLsizei, const GLuint*);
using PFNGLACTIVETEXTUREPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum);
using PFNGLBINDTEXTUREPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLuint);
using PFNGLTEXIMAGE2DPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLint, GLint, GLsizei, GLsizei,
                                                       GLint, GLenum, GLenum, const void*);
using PFNGLTEXPARAMETERIPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLenum, GLint);
using PFNGLGENFRAMEBUFFERSPROC = void (CHEMCAD_GL_APIENTRY*)(GLsizei, GLuint*);
using PFNGLDELETEFRAMEBUFFERSPROC = void (CHEMCAD_GL_APIENTRY*)(GLsizei, const GLuint*);
using PFNGLBINDFRAMEBUFFERPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLuint);
using PFNGLFRAMEBUFFERTEXTURE2DPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLenum, GLenum, GLenum, GLuint, GLint);
using PFNGLCHECKFRAMEBUFFERSTATUSPROC = GLenum (CHEMCAD_GL_APIENTRY*)(GLenum);
using PFNGLDRAWBUFFERSPROC = void (CHEMCAD_GL_APIENTRY*)(GLsizei, const GLenum*);
using PFNGLGENRENDERBUFFERSPROC = void (CHEMCAD_GL_APIENTRY*)(GLsizei, GLuint*);
using PFNGLDELETERENDERBUFFERSPROC = void (CHEMCAD_GL_APIENTRY*)(GLsizei, const GLuint*);
using PFNGLBINDRENDERBUFFERPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLuint);
using PFNGLRENDERBUFFERSTORAGEPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLenum, GLenum, GLsizei, GLsizei);
using PFNGLFRAMEBUFFERRENDERBUFFERPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLenum, GLenum, GLenum, GLuint);
using PFNGLDRAWARRAYSPROC = void (CHEMCAD_GL_APIENTRY*)(GLenum, GLint, GLsizei);
using PFNGLDRAWARRAYSINSTANCEDPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLenum, GLint, GLsizei, GLsizei);
using PFNGLDRAWELEMENTSPROC =
    void (CHEMCAD_GL_APIENTRY*)(GLenum, GLsizei, GLenum, const void*);

extern PFNGLGETSTRINGPROC glGetString;
extern PFNGLGETERRORPROC glGetError;
extern PFNGLGETINTEGERVPROC glGetIntegerv;
extern PFNGLGETBOOLEANVPROC glGetBooleanv;
extern PFNGLGETFLOATVPROC glGetFloatv;
extern PFNGLISENABLEDPROC glIsEnabled;
extern PFNGLENABLEPROC glEnable;
extern PFNGLDISABLEPROC glDisable;
extern PFNGLVIEWPORTPROC glViewport;
extern PFNGLCLEARPROC glClear;
extern PFNGLCLEARBUFFERFVPROC glClearBufferfv;
extern PFNGLDEPTHFUNCPROC glDepthFunc;
extern PFNGLDEPTHMASKPROC glDepthMask;
extern PFNGLBLENDFUNCSEPARATEPROC glBlendFuncSeparate;
extern PFNGLBLENDEQUATIONSEPARATEPROC glBlendEquationSeparate;
extern PFNGLCULLFACEPROC glCullFace;
extern PFNGLFRONTFACEPROC glFrontFace;
extern PFNGLGENBUFFERSPROC glGenBuffers;
extern PFNGLDELETEBUFFERSPROC glDeleteBuffers;
extern PFNGLBINDBUFFERPROC glBindBuffer;
extern PFNGLBUFFERDATAPROC glBufferData;
extern PFNGLGENVERTEXARRAYSPROC glGenVertexArrays;
extern PFNGLDELETEVERTEXARRAYSPROC glDeleteVertexArrays;
extern PFNGLBINDVERTEXARRAYPROC glBindVertexArray;
extern PFNGLENABLEVERTEXATTRIBARRAYPROC glEnableVertexAttribArray;
extern PFNGLCOLORMASKPROC glColorMask;
extern PFNGLCLEARCOLORPROC glClearColor;
extern PFNGLVERTEXATTRIBPOINTERPROC glVertexAttribPointer;
extern PFNGLVERTEXATTRIBDIVISORPROC glVertexAttribDivisor;
extern PFNGLCREATESHADERPROC glCreateShader;
extern PFNGLSHADERSOURCEPROC glShaderSource;
extern PFNGLCOMPILESHADERPROC glCompileShader;
extern PFNGLGETSHADERIVPROC glGetShaderiv;
extern PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
extern PFNGLDELETESHADERPROC glDeleteShader;
extern PFNGLCREATEPROGRAMPROC glCreateProgram;
extern PFNGLATTACHSHADERPROC glAttachShader;
extern PFNGLLINKPROGRAMPROC glLinkProgram;
extern PFNGLGETPROGRAMIVPROC glGetProgramiv;
extern PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog;
extern PFNGLDELETEPROGRAMPROC glDeleteProgram;
extern PFNGLUSEPROGRAMPROC glUseProgram;
extern PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
extern PFNGLUNIFORM1IPROC glUniform1i;
extern PFNGLUNIFORM1FPROC glUniform1f;
extern PFNGLUNIFORM2FPROC glUniform2f;
extern PFNGLUNIFORM3FPROC glUniform3f;
extern PFNGLUNIFORM4FPROC glUniform4f;
extern PFNGLUNIFORM1FVPROC glUniform1fv;
extern PFNGLUNIFORMMATRIX3FVPROC glUniformMatrix3fv;
extern PFNGLUNIFORMMATRIX4FVPROC glUniformMatrix4fv;
extern PFNGLGENTEXTURESPROC glGenTextures;
extern PFNGLDELETETEXTURESPROC glDeleteTextures;
extern PFNGLACTIVETEXTUREPROC glActiveTexture;
extern PFNGLBINDTEXTUREPROC glBindTexture;
extern PFNGLTEXIMAGE2DPROC glTexImage2D;
extern PFNGLTEXPARAMETERIPROC glTexParameteri;
extern PFNGLGENFRAMEBUFFERSPROC glGenFramebuffers;
extern PFNGLDELETEFRAMEBUFFERSPROC glDeleteFramebuffers;
extern PFNGLBINDFRAMEBUFFERPROC glBindFramebuffer;
extern PFNGLFRAMEBUFFERTEXTURE2DPROC glFramebufferTexture2D;
extern PFNGLCHECKFRAMEBUFFERSTATUSPROC glCheckFramebufferStatus;
extern PFNGLDRAWBUFFERSPROC glDrawBuffers;
extern PFNGLGENRENDERBUFFERSPROC glGenRenderbuffers;
extern PFNGLDELETERENDERBUFFERSPROC glDeleteRenderbuffers;
extern PFNGLBINDRENDERBUFFERPROC glBindRenderbuffer;
extern PFNGLRENDERBUFFERSTORAGEPROC glRenderbufferStorage;
extern PFNGLFRAMEBUFFERRENDERBUFFERPROC glFramebufferRenderbuffer;
extern PFNGLDRAWARRAYSPROC glDrawArrays;
extern PFNGLDRAWARRAYSINSTANCEDPROC glDrawArraysInstanced;
extern PFNGLDRAWELEMENTSPROC glDrawElements;

#undef CHEMCAD_GL_APIENTRY


using GlProcLoader = void* (*)(const char* name);

// Resolves every function below. Returns false if any required entry point is
// missing, in which case the renderer must stay disabled.
bool loadGl(GlProcLoader loader);
bool glLoaded();

// Human-readable driver identification, for the UI's diagnostics line.
const char* glVersionString();
const char* glRendererString();

}  // namespace chemcad::gfx
