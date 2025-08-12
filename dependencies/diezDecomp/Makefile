.SUFFIXES:
MAKEFLAGS += --no-builtin-rules --no-builtin-variables
SHELL=/bin/bash

FC := mpifort
FFLAGS :=
AR := ar rcs
LD := $(FC)
RM := rm -rf
GD := $(SRCS_DIR)/.gen-deps.awk
CPP := -cpp

ROOT_DIR := ../../
CONFIG_DIR := $(ROOT_DIR)/configs
BUILD_CONFIG_FILE ?= build.conf
include $(ROOT_DIR)/$(BUILD_CONFIG_FILE)
include $(CONFIG_DIR)/compilers.mk
include $(CONFIG_DIR)/flags.mk
include $(CONFIG_DIR)/libs.mk
ifeq ($(strip $(SINGLE_PRECISION)),1)
DEFINES += -D_DIEZDECOMP_SINGLE
endif

BUILD_DIR := ./build
SRC_DIR   := ./diezDecomp-repo/src
OBJ_DIR   := $(BUILD_DIR)
INC_DIR   := $(BUILD_DIR)/include
LIB_DIR   := $(BUILD_DIR)/lib

NAME := diezdecomp
LIB := $(patsubst %, lib%.a, $(NAME))
LIB := $(LIB_DIR)/$(LIB)
LIB  := $(LIB_DIR)/libdiezdecomp.a

SRC := $(wildcard $(SRC_DIR)/*.f90)

OBJ := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC))

all: $(LIB)

$(LIB): $(OBJ) | $(LIB_DIR) $(INC_DIR)
	$(AR) $@ $^
	cp $(BUILD_DIR)/*.mod $(INC_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) $(FFLAGS_MOD_DIR)$(BUILD_DIR) $(CPP) $(DEFINES) -c $< -o $@

$(OBJ_DIR)/diezdecomp_api_cans.o $(OBJ_DIR)/diezdecomp_api_generic.o: $(OBJ_DIR)/diezdecomp_core.o

$(OBJ_DIR):
	mkdir -p $@

$(LIB_DIR):
	mkdir -p $@

$(INC_DIR):
	mkdir -p $@

clean:
	$(RM) $(BUILD_DIR) *.o *.mod *.i

.PHONY: all clean
