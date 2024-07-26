#pragma once

// #include "../structure.h"

// namespace lible
// {
//     namespace ints
//     {
//         namespace GPUShells
//         {
//             struct ShellsData
//             {
//                 /*
//                  * Flattened data of shells to be used by GPU-kernels when calculating integrals.
//                  * Requires dimensions and offsets from 'ShellsInfo' to operate on.
//                  */
//                 double *contr_exps;
//                 double *contr_coeffs;
//                 double *norm_coeffs;
//                 int *cart_exps;
//             };

//             struct ShellsInfo
//             {
//                 /*
//                  * Information of each shell given by its index. Dimensions and offsets for
//                  * ShellsData when used by GPU integral kernels.
//                  */
//                 int *angmom;
//                 int *atomic_nr;
//                 int *dim_cart_exps;
//                 int *dim_contraction;
//                 int *dim_normalization;
//                 int *offset_cartexps;
//                 int *offset_contraction;
//                 int *offset_normalization;
//                 int *n_cartesian;
//                 int *n_spherical;
//                 int *pos_ao;
//             };
//         }

//         struct StructureGPU
//         {
//             /*
//              * Struct for representing the data necessary for calculating molecular integrals
//              * on the GPU.
//              */

//             StructureGPU(const std::unique_ptr<Structure> &structure);
//             ~StructureGPU();

//             double *allocateOneElInts(const size_t &n_ao);
//             void deallocateOneElInts(double *one_el_ints);
//             void copyOneElInts(const size_t &n_ao, double *one_el_ints_gpu, double *one_el_ints_cpu);

//             GPUShells::ShellsData *getShellsDataPtr()
//             {
//                 return &shells_data;
//             }

//             GPUShells::ShellsInfo *getShellsInfoPtr()
//             {
//                 return &shells_info;
//             }

//             std::vector<int> getShellCounts()
//             {
//                 return shell_counts;
//             }

//             std::vector<int> getShellOffsets()
//             {
//                 return shell_offsets;
//             }

//         private:
//             std::vector<int> shell_counts;
//             std::vector<int> shell_offsets;

//             GPUShells::ShellsData shells_data;
//             GPUShells::ShellsInfo shells_info;
//         };
//     }
// }