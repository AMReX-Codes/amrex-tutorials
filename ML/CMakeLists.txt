# LibTorch tutorials: -DAMReXTutorials_TORCH=ON, with Torch_DIR or CMAKE_PREFIX_PATH
# pointing to LibTorch (e.g., python3 -c "import torch; print(torch.utils.cmake_prefix_path)").
if(NOT AMReXTutorials_TORCH)
    message(STATUS "Skipping ML tutorials: AMReXTutorials_TORCH is OFF")
    return()
endif()

if(AMReX_GPU_BACKEND STREQUAL "HIP" OR AMReX_GPU_BACKEND STREQUAL "SYCL")
    message(STATUS "Skipping ML tutorials: only CPU and CUDA builds are supported")
    return()
endif()

find_package(Torch REQUIRED)

add_subdirectory(PyTorchInference)
if(TARGET ML_PyTorchInference)
    target_link_libraries(ML_PyTorchInference PRIVATE ${TORCH_LIBRARIES})
endif()
