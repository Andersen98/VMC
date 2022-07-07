env_base = Environment(CPPPATH=["/usr/include/eigen3", "./libigl/include", "./utils"])
serial = env_base.Clone(CCFLAGS=["-DSERIAL"])
serial.Library(
    "MicroVMC",
    ["CPSSlater.cpp"],
)
