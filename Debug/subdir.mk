################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../test_mathstuff.cpp \
../testmatrix.cpp \
../testvector.cpp 

OBJS += \
./test_mathstuff.o \
./testmatrix.o \
./testvector.o 

CPP_DEPS += \
./test_mathstuff.d \
./testmatrix.d \
./testvector.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


