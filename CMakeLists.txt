cmake_minimum_required(VERSION 3.14)
project(MethylSomaticAnalysis VERSION 1.0.0 LANGUAGES CXX)

# --------------------------------------------------
# 使用 FetchContent 自動取得 cxxopts (命令列參數解析庫)
include(FetchContent)
FetchContent_Declare(
  cxxopts
  GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
  GIT_TAG        v3.2.1       # 可依需求指定版本
)
FetchContent_MakeAvailable(cxxopts)
# --------------------------------------------------

# 設定 C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 編譯選項
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
  set(CMAKE_CXX_FLAGS_RELEASE   "${CMAKE_CXX_FLAGS_RELEASE}   -g -O3")
  set(CMAKE_CXX_FLAGS_DEBUG     "${CMAKE_CXX_FLAGS_DEBUG}     -g -O0 -DDEBUG")
endif()

# 輸出目錄
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# 依賴套件
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# 修改 htslib 檢測部分
find_path(HTSLIB_INCLUDE_DIR htslib/sam.h
  HINTS 
  $ENV{HTSLIB_ROOT}/include 
  /usr/include 
  /usr/local/include
)

find_library(HTSLIB_LIBRARY NAMES hts
  HINTS 
  $ENV{HTSLIB_ROOT}/lib 
  /usr/lib 
  /usr/local/lib
)

# 檢查 htslib 版本
if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
  # 設定內容檢查
  file(STRINGS "${HTSLIB_INCLUDE_DIR}/htslib/hts.h" HTSLIB_VERSION_DEFINE
    REGEX "#define HTS_VERSION .*$"
  )
  if(HTSLIB_VERSION_DEFINE MATCHES "#define HTS_VERSION \"([0-9]+\\.[0-9]+).*\"")
    set(HTSLIB_VERSION ${CMAKE_MATCH_1})
    message(STATUS "Found htslib version: ${HTSLIB_VERSION}")
    if(HTSLIB_VERSION VERSION_LESS "1.17")
      message(WARNING "htslib version >= 1.17 recommended (found ${HTSLIB_VERSION})")
    endif()
  else()
    message(WARNING "Could not determine htslib version")
  endif()
else()
  message(FATAL_ERROR "Could not find htslib")
endif()

# 如果 htslib < 1.17，添加自定義定義
if(HTSLIB_VERSION VERSION_LESS "1.17")
  add_definitions(-DHTSLIB_OLD_API)
endif()

# 包含目錄
include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${HTSLIB_INCLUDE_DIR}
  # cxxopts.hpp 已由 FetchContent 自動加入 include path
)

# 收集所有源碼
file(GLOB_RECURSE SOURCES "src/*.cpp")

# 建立可執行檔
add_executable(msa ${SOURCES})

# 連結所需函式庫
target_link_libraries(msa
  PRIVATE
    cxxopts::cxxopts       # cxxopts：命令列參數解析
    ${HTSLIB_LIBRARY}      # htslib：BAM/VCF 處理
    ZLIB::ZLIB             # zlib 壓縮
    Threads::Threads       # 多執行緒支援
)

# 安裝規則
install(TARGETS msa DESTINATION bin)

# 測試選項
option(BUILD_TESTS "Build tests" OFF)
if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# 配置摘要
message(STATUS "Build type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler:      ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "htslib include:    ${HTSLIB_INCLUDE_DIR}")
message(STATUS "htslib library:    ${HTSLIB_LIBRARY}")
message(STATUS "cxxopts version:   v3.2.1 (via FetchContent)")
message(STATUS "Build tests:       ${BUILD_TESTS}")
message(STATUS "Install prefix:    ${CMAKE_INSTALL_PREFIX}")
