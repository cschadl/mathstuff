################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../test_mathstuff.cpp \
../testgeom.cpp \
../testmatrix.cpp \
../testquaternion.cpp \
../testvector.cpp 

OBJS += \
./test_mathstuff.o \
./testgeom.o \
./testmatrix.o \
./testquaternion.o \
./testvector.o 

CPP_DEPS += \
./test_mathstuff.d \
./testgeom.d \
./testmatrix.d \
./testquaternion.d \
./testvector.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++11 -I/usr/include/eigen3 -O0 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


