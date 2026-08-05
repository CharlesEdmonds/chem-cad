#pragma once

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace chemcad::ui {

enum class FileDialogMode { Open, Save };

class FileDialog {
 public:
  void open(std::string title, FileDialogMode mode, std::vector<std::string> extensions,
            std::string defaultExtension = {}, std::string initialName = {});
  std::optional<std::string> draw();

 private:
  struct Entry {
    std::filesystem::path path;
    std::string label;
    bool directory = false;
  };

  void setDirectory(const std::filesystem::path& path);
  void refreshEntries();
  bool extensionAllowed(const std::filesystem::path& path) const;
  std::optional<std::string> accept();

  std::string title_;
  FileDialogMode mode_ = FileDialogMode::Open;
  std::vector<std::string> extensions_;
  std::string defaultExtension_;
  std::filesystem::path directory_;
  std::vector<Entry> entries_;
  std::array<char, 4096> pathBuffer_{};
  std::array<char, 512> filenameBuffer_{};
  std::string error_;
  bool open_ = false;
};

}  // namespace chemcad::ui
