################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../back/Copy\ of\ datas.cpp \
../back/Copy\ of\ main_kg.cpp 

OBJS += \
./back/Copy\ of\ datas.o \
./back/Copy\ of\ main_kg.o 

CPP_DEPS += \
./back/Copy\ of\ datas.d \
./back/Copy\ of\ main_kg.d 


# Each subdirectory must supply rules for building sources it contributes
back/Copy\ of\ datas.o: ../back/Copy\ of\ datas.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/opt/local/include/ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"back/Copy of datas.d" -MT"back/Copy\ of\ datas.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

back/Copy\ of\ main_kg.o: ../back/Copy\ of\ main_kg.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/opt/local/include/ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"back/Copy of main_kg.d" -MT"back/Copy\ of\ main_kg.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


