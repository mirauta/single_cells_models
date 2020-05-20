################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../additional_functions.cpp \
../datas.cpp \
../main_kg.cpp 

O_SRCS += \
../additional_functions.o \
../datas.o 

OBJS += \
./additional_functions.o \
./datas.o \
./main_kg.o 

CPP_DEPS += \
./additional_functions.d \
./datas.d \
./main_kg.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/opt/local/include/ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


