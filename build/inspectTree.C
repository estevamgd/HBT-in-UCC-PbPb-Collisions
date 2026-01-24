#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TKey.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TDirectory.h"

// Forward declaration for the recursive search function.
TTree* findTreeRecursive(TDirectory *dir);


/**
 * @brief Inspects a .root file to find the first TTree and lists its branches and their types.
 *
 * This macro opens the specified .root file and recursively searches through all
 * directories and subdirectories to find the first TTree. Once found, it prints the
 * tree's name and then lists all its branches (variables) and their data types.
 *
 * @usage
 * 1. Save this code as a file, for example, `inspectTree.C`.
 * 2. Open your terminal and run it with ROOT using the following command:
 * root -l 'inspectTree.C("your_file_name.root")'
 * (Replace "your_file_name.root" with the path to your actual file).
 *
 * @param filename The path to the .root file to inspect.
 */
void inspectTree(const char* filename) {
    // --- Open the ROOT file ---
    TFile *file = TFile::Open(filename, "READ");

    // --- File Opening Validation ---
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file '" << filename << "'" << std::endl;
        return;
    }

    std::cout << "Successfully opened file: " << filename << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    // --- Find the TTree using a recursive search ---
    TTree *tree = findTreeRecursive(file);

    // --- Tree Validation and Processing ---
    if (!tree) {
        std::cerr << "Error: No TTree found in file '" << filename << "'" << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Print the name of the TTree we found.
    std::cout << "Tree Name: " << tree->GetName() << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Variables (Branches) and Their Types:" << std::endl;
    std::cout << "----------------------------------------" << std::endl;


    // --- Get and List Branches ---
    TObjArray *branches = tree->GetListOfBranches();
    if (!branches) {
        std::cerr << "Error: Could not get list of branches." << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Iterate over all the branches.
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch *branch = dynamic_cast<TBranch*>(branches->At(i));
        if (branch) {
            TObjArray* leaves = branch->GetListOfLeaves();
            if (leaves && leaves->GetEntries() > 0) {
                TLeaf* leaf = (TLeaf*)leaves->At(0);
                if (leaf) {
                    printf("%-30s | Type: %s\n", branch->GetName(), leaf->GetTypeName());
                }
            } else {
                 printf("%-30s | Type: %s\n", branch->GetName(), branch->GetClassName());
            }
        }
    }
    std::cout << "----------------------------------------" << std::endl;

    // --- Cleanup ---
    file->Close();
    delete file;
}


/**
 * @brief Recursively searches a TDirectory for the first TTree.
 *
 * This function iterates through all keys in a given directory. If a key
 * corresponds to a TTree, it returns a pointer to it. If it corresponds
 * to another TDirectory, it calls itself on that subdirectory.
 *
 * @param dir The directory to search in.
 * @return A pointer to the first TTree found, or nullptr if none is found.
 */
TTree* findTreeRecursive(TDirectory *dir) {
    TIter next(dir->GetListOfKeys());
    TKey *key;
    TObject *obj;
    TTree *tree = nullptr;

    while ((key = (TKey*)next())) {
        obj = key->ReadObj();

        // --- ADDED NULL CHECK ---
        // If the object is unreadable (corrupted key), ReadObj() can return a nullptr.
        // Skip this key and continue to the next one to avoid a crash.
        if (!obj) {
            std::cerr << "Warning: Skipping unreadable object (key: " << key->GetName() << ")" << std::endl;
            continue;
        }

        // Check if the object is a TTree
        if (obj->InheritsFrom(TTree::Class())) {
            return dynamic_cast<TTree*>(obj);
        }
        // Check if the object is a TDirectory, and if so, search it
        else if (obj->InheritsFrom(TDirectory::Class())) {
            TDirectory *subdir = dynamic_cast<TDirectory*>(obj);
            tree = findTreeRecursive(subdir);
            // If we found a tree in the subdirectory, stop searching and return it.
            if (tree) {
                return tree;
            }
        }
    }
    // If we looped through all keys and found no tree, return nullptr.
    return nullptr;
}
