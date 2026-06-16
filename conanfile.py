from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout


class SteinhardtOrderParametersConan(ConanFile):
    name = "steinhardt-order-parameters"
    version = "1.0.4"
    package_type = "static-library"
    license = "MIT"
    settings = "os", "arch", "compiler", "build_type"
    requires = (
        "boost/1.88.0",
        "onetbb/2021.12.0",
        "coretoolkit/[>=2.5]",
        "spdlog/1.14.1",
        "nlohmann_json/3.11.3",
    )
    default_options = {"hwloc/*:shared": True}
    exports_sources = "CMakeLists.txt", "include/*", "src/*"

    def layout(self):
        cmake_layout(self)

    def generate(self):
        CMakeToolchain(self).generate()
        CMakeDeps(self).generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property(
            "cmake_target_name", "steinhardt-order-parameters::steinhardt-order-parameters"
        )
        self.cpp_info.libs = ["steinhardt-order-parameters_lib"]
        self.cpp_info.requires = [
            "boost::headers",
            "onetbb::onetbb",
            "coretoolkit::coretoolkit",
            "nlohmann_json::nlohmann_json",
            "spdlog::spdlog",
        ]
