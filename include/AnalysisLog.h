// AnalysisLog.h
#pragma once

#include <string>
#include <vector>
#include <ctime>
#include <fstream>
#include <filesystem>
#include <iomanip>

enum class LogEntryType {
    InputFile,
    SavedHistogram,
    SavedCanvas,
    SavedBenchmark,
    Generic
};

struct LogEntry {
    LogEntryType type;
    std::string description;
    std::string path;
    std::string filename;
    std::time_t timestamp;
};

// AnalysisLog.h (continued)

class AnalysisLog {
public:
    static AnalysisLog& instance() {
        static AnalysisLog logger;
        return logger;
    }

    // -------- INPUTS --------
    void addInputFile(const std::string& filepath) {
        entries.push_back({
            LogEntryType::InputFile,
            "Input ROOT file loaded",
            std::filesystem::path(filepath).parent_path().string(),
            std::filesystem::path(filepath).filename().string(),
            std::time(nullptr)
        });
    }

    // -------- OUTPUTS --------
    void addSavedObject(LogEntryType type,
                        const std::string& description,
                        const std::string& path,
                        const std::string& filename) {
        entries.push_back({
            type,
            description,
            path,
            filename,
            std::time(nullptr)
        });
    }

    // -------- FINAL SAVE --------
    void save(const std::string& path,
              const std::string& prefix) const {
        std::filesystem::create_directories(path);

        std::time_t t = std::time(nullptr);
        tm date = *std::localtime(&t);

        char fileName[256];
        sprintf(fileName,
                "%s/%s_log-%04d-%02d-%02d-%02d-%02d-%02d.txt",
                path.c_str(), prefix.c_str(),
                date.tm_year + 1900,
                date.tm_mon + 1,
                date.tm_mday,
                date.tm_hour,
                date.tm_min,
                date.tm_sec);

        std::ofstream out(fileName);
        out << "===== ANALYSIS LOG =====\n\n";

        for (const auto& e : entries) {
            tm ts = *std::localtime(&e.timestamp);

            out << "[" << std::put_time(&ts, "%Y-%m-%d %H:%M:%S") << "] ";

            switch (e.type) {
                case LogEntryType::InputFile:
                    out << "INPUT       ";
                    break;
                case LogEntryType::SavedHistogram:
                    out << "HISTOGRAM   ";
                    break;
                case LogEntryType::SavedCanvas:
                    out << "CANVAS      ";
                    break;
                case LogEntryType::SavedBenchmark:
                    out << "BENCHMARK   ";
                    break;
                default:
                    out << "INFO        ";
            }

            out << e.description << "\n";
            out << "    Path: " << e.path << "\n";
            out << "    File: " << e.filename << "\n\n";
        }

        out << "===== END OF LOG =====\n";
        out.close();
    }

private:
    AnalysisLog() = default;
    std::vector<LogEntry> entries;
};
