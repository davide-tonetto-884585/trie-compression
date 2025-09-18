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
#include <thread>
#include <mutex>
#include <future>
#include "../../include/XBWT/XBWT.hpp"
#include "../../include/ChainsDivisionSolver/ChainsDivisionSolver.hpp"
#include "../../include/TreeCompressor/TreeCompressor.hpp"

// Function to process a single tree file with given p value
struct CompressionResult {
    std::string filename;
    uint16_t p_value;
    uint64_t num_old_nodes;
    uint64_t num_new_nodes;
    uint64_t num_old_transitions;
    uint64_t num_new_transitions;
    uint64_t max_eq_classes;
};

// Structure to hold preprocessed tree data
struct PreprocessedTree {
    std::string filename;
    TreeDAWG<char> tree_dawg;
    std::vector<uint64_t> node_order;
    uint64_t num_old_nodes;
    uint64_t num_old_transitions;
    uint64_t max_eq_classes;
    
    PreprocessedTree(const std::string& filepath) {
        filename = std::filesystem::path(filepath).filename().string();
        
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

CompressionResult process_with_p_value(const PreprocessedTree& preprocessed, uint16_t p, bool verbose = false) {
    CompressionResult result;
    result.filename = preprocessed.filename;
    result.p_value = p;
    
    if (verbose)
        std::cout << "Processing " << result.filename << " with p=" << p << std::endl;
    
    // Build bipartite graph and solve (only operation that depends on p)
    ChainsDivisionSolver solver(preprocessed.tree_dawg, preprocessed.node_order, p, true, false);
    auto chains = solver.solve(false);
    
    // Validate that sum of chain lengths equals number of nodes
    uint64_t total_chain_length = 0;
    for (const auto& chain : chains) {
        total_chain_length += chain.size();
    }
    if (total_chain_length != preprocessed.num_old_nodes) {
        throw std::runtime_error("Chain validation failed: sum of chain lengths (" + 
                               std::to_string(total_chain_length) + 
                               ") != number of nodes (" + 
                               std::to_string(preprocessed.num_old_nodes) + ")");
    }
    
    // Compress tree
    TreeCompressor<char> compressor(preprocessed.tree_dawg, chains, false);
    auto transitions = compressor.get_transitions();
    
    // Fill result
    result.num_old_nodes = preprocessed.num_old_nodes;
    result.num_new_nodes = compressor.get_num_new_nodes();
    result.num_old_transitions = preprocessed.num_old_transitions;
    result.num_new_transitions = transitions.size();
    result.max_eq_classes = preprocessed.max_eq_classes;
    
    return result;
}

int main()
{
    bool verbose = false;
    std::string trees_directory = "./tree_generator/generated_trees/low";
    std::vector<uint16_t> p_values = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    
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
    std::cout << "Processing with p values: ";
    for (auto p : p_values) {
        std::cout << p << " ";
    }
    std::cout << std::endl;
    
    // Extract folder name from trees_directory for CSV filename
    std::string folder_name = std::filesystem::path(trees_directory).filename().string();
    std::string csv_filename = "compression_results_" + folder_name + ".csv";
    
    // Open CSV file for writing results
    std::ofstream csv_file(csv_filename);
    csv_file << "filename,p_value,num_old_nodes,num_new_nodes,num_old_transitions,num_new_transitions,max_eq_classes\n";
    
    // Process each tree file with each p value
    int total_experiments = tree_files.size() * p_values.size();
    int current_experiment = 0;
    std::mutex csv_mutex; // Mutex for thread-safe CSV writing
    std::mutex progress_mutex; // Mutex for thread-safe progress updates
    
    for (const auto& tree_file : tree_files) {
        try {
            // Preprocess tree once (expensive operations)
            std::cout << "Preprocessing " << std::filesystem::path(tree_file).filename().string() << "..." << std::endl;
            PreprocessedTree preprocessed(tree_file);
            
            // Create futures for parallel processing of p values
            std::vector<std::future<CompressionResult>> futures;
            
            // Launch threads for each p value
            for (auto p : p_values) {
                futures.push_back(std::async(std::launch::async, [&preprocessed, p, verbose]() {
                    return process_with_p_value(preprocessed, p, verbose);
                }));
            }
            
            // Collect results from all threads
            for (size_t i = 0; i < futures.size(); ++i) {
                try {
                    CompressionResult result = futures[i].get();
                    
                    // Thread-safe progress update
                    {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        current_experiment++;
                        std::cout << "[" << current_experiment << "/" << total_experiments << "] ";
                        std::cout << "Completed " << result.filename << " (p=" << result.p_value << ")" << std::endl;
                    }
                    
                    // Thread-safe CSV writing
                    {
                        std::lock_guard<std::mutex> lock(csv_mutex);
                        csv_file << result.filename << ","
                                << result.p_value << ","
                                << result.num_old_nodes << ","
                                << result.num_new_nodes << ","
                                << result.num_old_transitions << ","
                                << result.num_new_transitions << ","
                                << result.max_eq_classes << "\n";
                    }
                }
                catch (const std::exception& e) {
                    std::lock_guard<std::mutex> lock(progress_mutex);
                    current_experiment++;
                    std::cerr << "Error processing with p=" << p_values[i] << ": " << e.what() << std::endl;
                }
            }
        }
        catch (const std::exception& e) {
            std::cerr << "Error preprocessing " << tree_file << ": " << e.what() << std::endl;
            // Skip all p values for this file if preprocessing fails
            current_experiment += p_values.size();
        }
    }
    
    csv_file.close();
    std::cout << "\nAll experiments completed! Results saved to compression_results.csv" << std::endl;
    
    return 0;
}
