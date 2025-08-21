Below is a brief guide to the contents of the `MIP_Feature_Computation` folder (often called “mip\_feature\_computation” in Hutter et al.’s repository). This directory’s purpose is to build and run the code that extracts numerical features from every MIP instance—features like numbers of variables/constraints, constraint‐matrix density, presolve statistics, and so on—that feed into their empirical performance models.

---

## Top‐Level Executables & Scripts

1. **`MIPfeature`**

   * This is the primary executable (a compiled binary) that you run to compute a full set of features for a given MIP instance. Under the hood, it links together various C++ “comp\_feat” files to collect statistics at presolve time, LP relaxation time, and simple graph‐based measures.
   * In practice, you invoke it on a file (e.g. `MIPfeature some_instance.mps`) and it prints out a line of comma‐separated features (matching the descriptions in `FeatureDescription.txt`).

2. **`MIPfeature‐used‐for‐almost‐everything`**

   * This is essentially a symlink or alternate name for the same “full” feature extractor. Both `MIPfeature` and this file point to the same compiled binary; the extra name reminds you that this is the version used in almost every feature computation pipeline (including low‐memory or specialized variants).

3. **`Makefile`**

   * Contains the build rules to compile all the `.cpp` files into `MIPfeature` (and other test binaries). Simply running `make` in this directory builds the main executable(s).

---

## Source Files for Feature Extraction

All of the files whose names begin with `comp_feat…` are different versions or experiments of the feature‐extraction code. Each one implements more or fewer features, or uses alternative techniques to reduce memory footprint. In particular:

1. **`comp_feat.cpp`**

   * The “canonical” C++ source that defines the bulk of Hutter et al.’s MIP feature set. When you compile the directory with the provided `Makefile`, this file (along with its headers) is the main driver behind `MIPfeature`.
   * It orchestrates calls to CPLEX (via their C API) to collect presolve statistics (e.g., number of rows/columns eliminated, number of cuts added) and then computes combinatorial features (like bipartite‐graph node degrees) once the presolve is finished.

2. **`comp_feat.h`**

   * Header file where all feature definitions, data‐structures, and function prototypes live. This is included by every `.cpp` that needs to reference the feature‐computation routines.

3. **`comp_feat_lin.cpp`**

   * A “linear‐only” version of the feature code: it collects only those features that can be computed from the **pure LP relaxation** (that is, after ignoring integrality). This version omits any MIP‐specific counts (no integer‐var classes, no cut rounds). It is faster if you only care about LP‐based features.

4. **`lin_comp_feat.cpp`**

   * Very similar to `comp_feat_lin.cpp`—another take on extracting LP‐specific features. In some historical checks, Hutter’s team experimented with two different LP‐extraction approaches; both produce feature sets restricted to linear‐program metrics (e.g., basis sparsity, objective‐function coefficients).

5. **`comp_feat_low_mem.cpp`**

   * A “low‐memory” variant that trades off some computed features in order to reduce peak memory usage. If you want to compute features on very large MIPs without running out of RAM, compile and use this version instead of the full `comp_feat.cpp`. It drops or simplifies any feature whose computation requires building large intermediate data structures (for instance, dense adjacency arrays).

6. **`comp_feat_try1.cpp`**

   * An early prototype (“trial 1”) of the feature extractor. It may compute just a subset of the final features, or use a different traversal order over CPLEX’s internal data structures. It mostly serves as a reference to show how the feature set evolved.

7. **`orig_comp_feat.cpp`**, **`old_comp_feat.cpp`**, **`new_orig.cpp`**, **`frank_save_comp_feat.cpp`**

   * These are older or alternative versions of the same feature‐computation code. Over time, Hutter et al. tried different implementations (for example, to stabilize numerical results or to save intermediate files in a particular format).
   * You generally only need to worry about `comp_feat.cpp` (and its low‐mem or linear variants) unless you’re curious about how the code was refactored over various paper revisions.

---

## Utility & Helper Files

1. **`set_example.cpp`**, **`set2_example.cpp`**

   * Small demonstration programs that show how to invoke `comp_feat.cpp`’s internals. If you want to write your own wrapper (for instance, to embed feature extraction in a Python script), these example files illustrate the minimal calls needed (e.g., reading an MPS file, calling `computeFeatures()`, printing results).

2. **`stopwatch.cpp`**, **`stopwatch.h`**, **`stopwatch.o`**

   * A simple timing utility used throughout the feature‐computation code. The feature extractor measures how long presolve took, how long LP solves took, etc., using these stopwatch routines.

3. **`cplex.log`**

   * A sample log file showing one run of CPLEX’s initialization and presolve. It’s mostly for debugging: if you want to see exactly which CPLEX calls were issued (e.g. `CPXpresolve()`, `CPXlpopt()`), you can refer to this log.

4. **`README.txt`**

   * Provides high‐level instructions on how to build and invoke the feature extractor. It lists required dependencies (e.g. CPLEX include/lib paths) and gives example compile/run commands.

5. **`FeatureDescription.txt`**

   * For each numerical column that `MIPfeature` outputs, this text file explains what the column means. For example,

     * *“nCols”* = total number of columns in the original MIP,
     * *“avgRowDensity”* = average number of nonzeros per constraint row,
     * *“PresolvedRows”* = number of rows remaining after CPLEX’s presolve.
   * It’s the definitive reference if you want to know exactly which feature index corresponds to which metric.

6. **`instancelist.txt`**

   * A plain‐text file containing a list of all MIP instance filenames (one per line) for which you typically compute features. Often, Hutter’s scripts iterate through this list, invoking `MIPfeature` on each `.mps` file.

---

## Miscellaneous Source Files

1. **`old_comp_feat.cpp`**, **`orig_comp_feat.cpp`**, **`new_orig.cpp`**

   * These are historical or experimental versions of the same feature code. They exist largely for archival purposes to show how the feature suite was refined. Typically you won’t compile or use them in a new workflow.

2. **`comp_feat.o`**, **`stopwatch.o`**

   * Object files generated by previous compiles. If you run `make clean` and then `make`, these get rebuilt. You can usually ignore `.o` files unless you want to inspect the build process.

---

## How to Use It in Practice

1. **Compile everything**

   ```bash
   cd MIP_Feature_Computation
   make
   ```

   This produces the `MIPfeature` (and related) binaries.

2. **Compute features on a single MIP**

   ```bash
   ./MIPfeature path/to/my_instance.mps > my_instance.features.csv
   ```

   The output will be a single CSV line:

   ```
   filename, nCols, nRows, avgRowDensity, PresolvedRows, ..., TimingStats...
   ```

   * If you want only LP‐based features, compile `comp_feat_lin.cpp` (change the `Makefile` to use it) and run that binary instead.

3. **Batch‐compute features for a list of instances**
   Suppose you have a file called `instancelist.txt` containing hundreds of `.mps` filenames. You can loop over them in a shell script:

   ```bash
   while read instance; do
     ./MIPfeature "$instance" >> all_instances.features.csv
   done < instancelist.txt
   ```

   This produces one big CSV with one row per instance.

4. **Use low‐memory mode**
   If your MIPs are huge and the standard `MIPfeature` runs out of RAM, open the `Makefile` and replace references to `comp_feat.cpp` with `comp_feat_low_mem.cpp`. Then `make` again. The resulting binary will skip or approximate some high‐cost features so that it can handle larger models.

5. **Consult `FeatureDescription.txt`**
   After you generate a CSV (say `all_instances.features.csv`), open `FeatureDescription.txt` to see what each column index means (e.g., column 5 = number of integer variables, column 12 = maximum variable degree, etc.). You need this legend when you later load the features into a machine‐learning pipeline.

---

### In Summary

* **`MIPfeature`** (and its alias `MIPfeature‐used‐for‐almost‐everything`) is the compiled program you actually run to extract features.
* All the `comp_feat*.cpp` files are different implementations or variants of that feature‐extraction logic (full‐memory vs low‐memory vs LP‐only).
* The header `comp_feat.h` and helper files (`stopwatch.*`, `FeatureDescription.txt`) define exactly which features are computed and how.
* Utility examples (`set_example.cpp`) show minimal code snippets for calling the main routines if you want to embed feature computation elsewhere.
* Finally, `instancelist.txt` is the canonical list of MIP files to process, and the `.log` file is for debugging CPLEX calls.

By compiling and running the code in this folder, you produce exactly the same per-instance feature vectors that Hutter et al. used to train their empirical performance models for CPLEX, Gurobi, SCIP, and other solvers.
