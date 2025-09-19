#include <iostream>
#include <TreeDAWG.hpp>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <chrono>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <thread>
#include <mutex>
#include <future>
#include "../../include/XBWT/XBWT.hpp"
#include "../../include/ChainsDivisionSolver/ChainsDivisionSolver.hpp"
#include "../../include/TreeCompressor/TreeCompressor.hpp"

// Structure to hold time measurement results
struct TimeResult {
    std::string filename;
    uint16_t p_value;
    uint64_t num_old_nodes;
    uint64_t num_new_nodes_optimized_false;
    uint64_t num_new_nodes_optimized_true;
    uint64_t num_old_transitions;
    uint64_t num_new_transitions_optimized_false;
    uint64_t num_new_transitions_optimized_true;
    uint64_t max_eq_classes;
    double time_optimized_false_ms;  // Time in milliseconds for optimized=false
    double time_optimized_true_ms;   // Time in milliseconds for optimized=true
};

// Structure to hold preprocessed tree data
struct PreprocessedTree {
    std::string filename;
    TreeDAWG<char> tree_dawg;
    std::vector<uint64_t> node_order;
    uint64_t num_old_nodes;
    uint64_t num_old_transitions;
    uint64_t max_eq_classes;
    
    PreprocessedTree(const std::string& filepath) : filename(std::filesystem::path(filepath).filename().string()) {
        
        // Read tree from file
        std::string str;
        std::ifstream inputFile(filepath);
        if (!inputFile.is_open()) {
            throw std::runtime_error("Unable to open file: " + filepath);
        }
        getline(inputFile, str);
        inputFile.close();
        
        // Parse tree
        LabeledTree<char> tree(str, [](const std::string &label) -> char { return label[0]; });
        
        // Build Tree DAWG (expensive operation - done once)
        tree_dawg = TreeDAWG<char>(tree);
        tree_dawg.minimize();
        
        // Compute node order (expensive operation - done once)
        node_order = tree_dawg.stable_sort_nodes_by_label_path(false);
        
        // Store frequently used values
        num_old_nodes = tree_dawg.get_num_nodes();
        num_old_transitions = tree_dawg.get_num_transitions();

        max_eq_classes = 0;
        for (uint64_t i = 0; i < num_old_nodes; ++i)
        {
            if (max_eq_classes < tree_dawg[i].get_equivalence_class())
                max_eq_classes = tree_dawg[i].get_equivalence_class();
        }
    }
};

// Function to process a single p_value for a given tree (thread-safe)
// Measures time 5 times and returns the average for more accurate results
TimeResult process_p_value_parallel(const PreprocessedTree& preprocessed, uint16_t p_value) {
    TimeResult result;
    
    const int num_runs = 5;
    std::vector<double> times_optimized_false;
    std::vector<double> times_optimized_true;
    
    // Run the experiment 5 times to get more accurate timing
    for (int run = 0; run < num_runs; ++run) {
        // Measure time for optimized=false
        {
            auto start_time = std::chrono::high_resolution_clock::now();
            ChainsDivisionSolver solver_false(preprocessed.tree_dawg, preprocessed.node_order, p_value, false, false);
            auto chains_false = solver_false.solve(false);
            
            // Compress tree to get final metrics for optimized=false (only on first run)
            if (run == 0) {
                TreeCompressor<char> compressor_false(preprocessed.tree_dawg, chains_false, false);
                auto transitions_false = compressor_false.get_transitions();
                result.num_new_nodes_optimized_false = compressor_false.get_num_new_nodes();
                result.num_new_transitions_optimized_false = transitions_false.size();
            }
            
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration_false = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            times_optimized_false.push_back(duration_false.count() / 1000.0);
        }
        
        // Measure time for optimized=true
        {
            auto start_time = std::chrono::high_resolution_clock::now();
            ChainsDivisionSolver solver_true(preprocessed.tree_dawg, preprocessed.node_order, p_value, true, false);
            auto chains_true = solver_true.solve(false);
            
            // Compress tree to get final metrics for optimized=true (only on first run)
            if (run == 0) {
                TreeCompressor<char> compressor_true(preprocessed.tree_dawg, chains_true, false);
                auto transitions_true = compressor_true.get_transitions();
                result.num_new_nodes_optimized_true = compressor_true.get_num_new_nodes();
                result.num_new_transitions_optimized_true = transitions_true.size();
            }
            
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration_true = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            times_optimized_true.push_back(duration_true.count() / 1000.0);
        }
    }
    
    // Calculate averages
    double avg_time_false = 0.0;
    double avg_time_true = 0.0;
    
    for (double time : times_optimized_false) {
        avg_time_false += time;
    }
    avg_time_false /= num_runs;
    
    for (double time : times_optimized_true) {
        avg_time_true += time;
    }
    avg_time_true /= num_runs;
    
    // Fill common result fields
    result.filename = preprocessed.filename;
    result.p_value = p_value;
    result.num_old_nodes = preprocessed.num_old_nodes;
    result.num_old_transitions = preprocessed.num_old_transitions;
    result.max_eq_classes = preprocessed.max_eq_classes;
    result.time_optimized_false_ms = avg_time_false;
    result.time_optimized_true_ms = avg_time_true;
    
    return result;
}

int main()
{
    bool verbose = true;
    std::string trees_directory = "./tree_generator/generated_trees/time_exp";
    std::vector<uint16_t> p_values = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    
    // Get all .txt files in the directory
    std::vector<std::string> tree_files;
    for (const auto& entry : std::filesystem::directory_iterator(trees_directory)) {
        if (entry.path().extension() == ".txt") {
            tree_files.push_back(entry.path().string());
        }
    }
    
    // Sort files for consistent ordering
    std::sort(tree_files.begin(), tree_files.end());
    
    std::cout << "Found " << tree_files.size() << " tree files to process" << std::endl;
    std::cout << "Testing p values: ";
    for (auto p : p_values) {
        std::cout << p << " ";
    }
    std::cout << std::endl;
    
    // Create CSV filename with timestamp
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t), "%Y%m%d_%H%M%S");
    std::string csv_filename = "time_experiments_" + ss.str() + ".csv";
    
    // Open CSV file for writing results
    std::ofstream csv_file(csv_filename);
    csv_file << "filename,p_value,num_old_nodes,num_new_nodes_optimized_false,num_new_nodes_optimized_true,num_old_transitions,num_new_transitions_optimized_false,num_new_transitions_optimized_true,max_eq_classes,time_optimized_false_ms,time_optimized_true_ms\n";
    
    // Process each tree file with parallel p-value computation
    int total_experiments = tree_files.size() * p_values.size();
    int current_experiment = 0;
    
    for (const auto& tree_file : tree_files) {
        try {
            // Preprocess tree once (expensive operations)
            std::cout << "\nPreprocessing " << std::filesystem::path(tree_file).filename().string() << "..." << std::endl;
            PreprocessedTree preprocessed(tree_file);
            
            std::cout << "Tree info: " << preprocessed.num_old_nodes << " nodes, " 
                      << preprocessed.num_old_transitions << " transitions, " 
                      << preprocessed.max_eq_classes + 1 << " equivalence classes" << std::endl;
            
            // Determine number of threads (use hardware concurrency or default to 4)
            unsigned int num_threads = std::thread::hardware_concurrency();
            if (num_threads == 0) num_threads = 4;
            
            std::cout << "  Using " << num_threads << " threads for p-values (5 runs each for accurate timing)" << std::endl;
            
            // Process p_values in parallel using futures
            std::vector<std::future<TimeResult>> futures;
            
            for (auto p : p_values) {
                futures.push_back(std::async(std::launch::async, 
                    process_p_value_parallel, std::cref(preprocessed), p));
            }
            
            // Collect results and maintain order
            std::vector<TimeResult> tree_results;
            for (size_t i = 0; i < futures.size(); ++i) {
                try {
                    current_experiment++;
                    
                    TimeResult result = futures[i].get();
                    tree_results.push_back(result);
                    
                    std::cout << "[" << current_experiment << "/" << total_experiments << "] "
                             << "Completed " << result.filename << " (p=" << result.p_value 
                             << ") - optimized=false: " << std::fixed << std::setprecision(2) 
                             << result.time_optimized_false_ms << "ms (avg of 5 runs), "
                             << "optimized=true: " << result.time_optimized_true_ms << "ms (avg of 5 runs)" << std::endl;
                }
                catch (const std::exception& e) {
                    std::cerr << "Error processing " << tree_file << " with p=" << p_values[i] << ": " << e.what() << std::endl;
                }
            }
            
            // Sort results by p_value to maintain consistent order in CSV
            std::sort(tree_results.begin(), tree_results.end(), 
                     [](const TimeResult& a, const TimeResult& b) {
                         return a.p_value < b.p_value;
                     });
            
            // Write all results for this tree to CSV
            for (const auto& result : tree_results) {
                csv_file << result.filename << ","
                        << result.p_value << ","
                        << result.num_old_nodes << ","
                        << result.num_new_nodes_optimized_false << ","
                        << result.num_new_nodes_optimized_true << ","
                        << result.num_old_transitions << ","
                        << result.num_new_transitions_optimized_false << ","
                        << result.num_new_transitions_optimized_true << ","
                        << result.max_eq_classes << ","
                        << std::fixed << std::setprecision(3) << result.time_optimized_false_ms << ","
                        << std::fixed << std::setprecision(3) << result.time_optimized_true_ms << "\n";
            }
            
            // Flush CSV file to avoid data loss
            csv_file.flush();
        }
        catch (const std::exception& e) {
            std::cerr << "Error preprocessing " << tree_file << ": " << e.what() << std::endl;
            // Skip all p values for this file if preprocessing fails
            current_experiment += p_values.size();
        }
    }
    
    csv_file.close();
    std::cout << "\nAll time experiments completed! Results saved to " << csv_filename << std::endl;
    
    return 0;
}