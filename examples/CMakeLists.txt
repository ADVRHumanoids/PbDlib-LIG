cmake_minimum_required(VERSION 2.8)
PROJECT(pbdlib)

add_executable(test_gmm test_gmm.cpp)
target_link_libraries(test_gmm pbd)

add_executable(test_gmr test_gmr.cpp)
target_link_libraries(test_gmr pbd)

add_executable(test_gmr2 test_gmr2.cpp)
target_link_libraries(test_gmr2 pbd)

add_executable(test_tpgmm_timeInvariant test_tpgmm_timeInvariant.cpp)
target_link_libraries(test_tpgmm_timeInvariant pbd)

add_executable(test_lqr test_lqr.cpp)
target_link_libraries(test_lqr pbd)

add_executable(test_onlineDP test_onlineDP.cpp)
target_link_libraries(test_onlineDP pbd)

add_executable(test_rewardWeightedRefinement test_rewardWeightedRefinement.cpp)
target_link_libraries(test_rewardWeightedRefinement pbd)

add_executable(test_datapoints test_datapoints.cpp)
target_link_libraries(test_datapoints pbd)

add_executable(test_demonstration test_demonstration.cpp)
target_link_libraries(test_demonstration pbd)

add_executable(test_tpdemonstration test_tpdemonstration.cpp)
target_link_libraries(test_tpdemonstration pbd)

add_executable(test_mvn test_mvn.cpp)
target_link_libraries(test_mvn pbd)

add_executable(test_HSMM test_hsmm.cpp)
target_link_libraries(test_HSMM pbd)

add_executable(test_HSMM2 test_hsmm2.cpp)
target_link_libraries(test_HSMM2 pbd)

add_executable(test_ADHSMM test_adhsmm.cpp)
target_link_libraries(test_ADHSMM pbd)

add_executable(test_quaternions test_quaternions.cpp)
target_link_libraries(test_quaternions pbd)

add_executable(test_tpgmm_quaternions test_tpgmm_quaternions.cpp)
target_link_libraries(test_tpgmm_quaternions pbd)

add_executable(test_TrajGMM test_trajGMM.cpp)
target_link_libraries(test_TrajGMM pbd)

add_executable(test_mpc test_mpc.cpp)
target_link_libraries(test_mpc pbd)

add_executable(test_trajMPC test_trajMPC.cpp)
target_link_libraries(test_trajMPC pbd)
