#include "ui/file_dialog.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>

#include "imgui.h"

namespace chemcad::ui {
namespace {

std::filesystem::path& lastDirectory() {
  static std::filesystem::path path;
  return path;
}

std::string lower(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

template <size_t N>
void setBuffer(std::array<char, N>& buffer, const std::string& text) {
  std::snprintf(buffer.data(), buffer.size(), "%s", text.c_str());
}

std::filesystem::path initialDirectory() {
  if (!lastDirectory().empty()) return lastDirectory();
  if (const char* home = std::getenv("HOME"); home && *home) return home;
  std::error_code ec;
  const auto cwd = std::filesystem::current_path(ec);
  return ec ? std::filesystem::path(".") : cwd;
}

}  // namespace

void FileDialog::open(std::string title, FileDialogMode mode,
                      std::vector<std::string> extensions, std::string defaultExtension,
                      std::string initialName) {
  title_ = std::move(title);
  mode_ = mode;
  extensions_.clear();
  extensions_.reserve(extensions.size());
  for (std::string ext : extensions) {
    ext = lower(std::move(ext));
    if (!ext.empty() && ext.front() != '.') ext.insert(ext.begin(), '.');
    extensions_.push_back(std::move(ext));
  }
  defaultExtension_ = lower(std::move(defaultExtension));
  if (!defaultExtension_.empty() && defaultExtension_.front() != '.') {
    defaultExtension_.insert(defaultExtension_.begin(), '.');
  }
  setBuffer(filenameBuffer_, initialName);
  error_.clear();
  open_ = true;
  setDirectory(initialDirectory());
  ImGui::OpenPopup(title_.c_str());
}

bool FileDialog::extensionAllowed(const std::filesystem::path& path) const {
  if (extensions_.empty()) return true;
  const std::string ext = lower(path.extension().string());
  return std::find(extensions_.begin(), extensions_.end(), ext) != extensions_.end();
}

void FileDialog::setDirectory(const std::filesystem::path& path) {
  std::error_code ec;
  auto candidate = std::filesystem::absolute(path, ec);
  if (ec) candidate = path;
  candidate = candidate.lexically_normal();
  if (!std::filesystem::is_directory(candidate, ec) || ec) {
    error_ = "Not a readable directory: " + candidate.string();
    return;
  }
  directory_ = candidate;
  lastDirectory() = candidate;
  setBuffer(pathBuffer_, directory_.string());
  refreshEntries();
}

void FileDialog::refreshEntries() {
  entries_.clear();
  error_.clear();

  if (directory_.has_parent_path() && directory_ != directory_.root_path()) {
    entries_.push_back({directory_.parent_path(), "../", true});
  }

  std::vector<Entry> children;
  std::error_code ec;
  std::filesystem::directory_iterator it(directory_,
                                         std::filesystem::directory_options::skip_permission_denied,
                                         ec);
  if (ec) {
    error_ = "Cannot read directory: " + ec.message();
    return;
  }
  const std::filesystem::directory_iterator end;
  while (it != end) {
    const auto path = it->path();
    std::error_code typeError;
    const bool directory = it->is_directory(typeError);
    if (!typeError && (directory || extensionAllowed(path))) {
      std::string label = path.filename().string();
      if (directory) label += '/';
      children.push_back({path, std::move(label), directory});
    }
    it.increment(ec);
    if (ec) {
      error_ = "Directory listing stopped: " + ec.message();
      break;
    }
  }

  std::sort(children.begin(), children.end(), [](const Entry& a, const Entry& b) {
    if (a.directory != b.directory) return a.directory > b.directory;
    return lower(a.label) < lower(b.label);
  });
  entries_.insert(entries_.end(), std::make_move_iterator(children.begin()),
                  std::make_move_iterator(children.end()));
}

std::optional<std::string> FileDialog::accept() {
  if (filenameBuffer_[0] == '\0') {
    error_ = "Choose a file name.";
    return std::nullopt;
  }

  std::filesystem::path path(filenameBuffer_.data());
  if (path.is_relative()) path = directory_ / path;
  if (mode_ == FileDialogMode::Save && path.extension().empty() && !defaultExtension_.empty()) {
    path += defaultExtension_;
  }
  path = path.lexically_normal();

  if (!extensionAllowed(path)) {
    error_ = "The selected file does not match this dialog's file type.";
    return std::nullopt;
  }

  std::error_code ec;
  if (mode_ == FileDialogMode::Open) {
    if (!std::filesystem::is_regular_file(path, ec) || ec) {
      error_ = "File does not exist or is not readable.";
      return std::nullopt;
    }
  } else {
    const auto parent = path.has_parent_path() ? path.parent_path() : directory_;
    if (!std::filesystem::is_directory(parent, ec) || ec) {
      error_ = "The destination directory does not exist.";
      return std::nullopt;
    }
  }

  if (path.has_parent_path()) lastDirectory() = path.parent_path();
  return path.string();
}

std::optional<std::string> FileDialog::draw() {
  if (!open_) return std::nullopt;

  ImGui::SetNextWindowSize(ImVec2(620.0f, 470.0f), ImGuiCond_Appearing);
  bool keepOpen = true;
  std::optional<std::string> result;
  if (ImGui::BeginPopupModal(title_.c_str(), &keepOpen, ImGuiWindowFlags_NoCollapse)) {
    ImGui::SetNextItemWidth(-1.0f);
    if (ImGui::InputText("Path", pathBuffer_.data(), pathBuffer_.size(),
                         ImGuiInputTextFlags_EnterReturnsTrue)) {
      setDirectory(pathBuffer_.data());
    }

    if (ImGui::BeginListBox("##directory", ImVec2(-1.0f, 300.0f))) {
      for (const Entry& entry : entries_) {
        const bool selected = !entry.directory && entry.path.filename() == filenameBuffer_.data();
        if (ImGui::Selectable(entry.label.c_str(), selected,
                              entry.directory ? ImGuiSelectableFlags_AllowDoubleClick : 0)) {
          if (entry.directory) {
            if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left) || entry.label == "../") {
              setDirectory(entry.path);
              break;
            }
          } else {
            setBuffer(filenameBuffer_, entry.path.filename().string());
          }
        }
      }
      ImGui::EndListBox();
    }

    ImGui::SetNextItemWidth(-1.0f);
    const bool enter = ImGui::InputText("File name", filenameBuffer_.data(), filenameBuffer_.size(),
                                        ImGuiInputTextFlags_EnterReturnsTrue);
    if (!error_.empty()) ImGui::TextColored(ImVec4(0.95f, 0.40f, 0.35f, 1.0f), "%s", error_.c_str());

    const char* action = mode_ == FileDialogMode::Open ? "Open" : "Save";
    if (enter || ImGui::Button(action)) {
      result = accept();
      if (result) {
        open_ = false;
        ImGui::CloseCurrentPopup();
      }
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel")) {
      open_ = false;
      ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
  }

  if (!keepOpen) open_ = false;
  return result;
}

}  // namespace chemcad::ui
