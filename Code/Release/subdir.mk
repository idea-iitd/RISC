################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../AOR.c \
../AOR_BF.c \
../arrayListInt.c \
../arrayListInt_Flexible.c \
../binaryFingerPrint.c \
../chemFP.c \
../data.c \
../dataBinary.c \
../dataNonBinary.c \
../divideSkip.c \
../divideSkip_BF.c \
../experiments.c \
../helper.c \
../invertedIndexSmall.c \
../invertedIndex_BF.c \
../minHeap.c \
../minHeapInt.c \
../minHeapPlain.c \
../nonBinaryFingerPrint.c \
../options.c \
../risc.c 

OBJS += \
./AOR.o \
./AOR_BF.o \
./arrayListInt.o \
./arrayListInt_Flexible.o \
./binaryFingerPrint.o \
./chemFP.o \
./data.o \
./dataBinary.o \
./dataNonBinary.o \
./divideSkip.o \
./divideSkip_BF.o \
./experiments.o \
./helper.o \
./invertedIndexSmall.o \
./invertedIndex_BF.o \
./minHeap.o \
./minHeapInt.o \
./minHeapPlain.o \
./nonBinaryFingerPrint.o \
./options.o \
./risc.o 

C_DEPS += \
./AOR.d \
./AOR_BF.d \
./arrayListInt.d \
./arrayListInt_Flexible.d \
./binaryFingerPrint.d \
./chemFP.d \
./data.d \
./dataBinary.d \
./dataNonBinary.d \
./divideSkip.d \
./divideSkip_BF.d \
./experiments.d \
./helper.d \
./invertedIndexSmall.d \
./invertedIndex_BF.d \
./minHeap.d \
./minHeapInt.d \
./minHeapPlain.d \
./nonBinaryFingerPrint.d \
./options.d \
./risc.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -std=c99 -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


