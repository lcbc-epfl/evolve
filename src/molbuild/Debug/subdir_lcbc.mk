################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../MoleculeBuilder.cpp \
../main.cpp 

OBJS += \
./MoleculeBuilder.o \
./main.o

CPP_DEPS += \
./MoleculeBuilder.d \
./main.d


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	clang++ -I/software/boost-1.55.0/include -I/software/openbabel-2.3.2/include/openbabel-2.0 -O3 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


