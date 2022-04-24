#
# adapted from: https://fortran-lang.org/learn/building_programs/project_make
#
# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Project name
NAME := cans

TARGET := $(NAME)

PWD=$(shell pwd)
ROOT_DIR := $(PWD)
SRC_DIR := $(ROOT_DIR)/src
APP_DIR := $(ROOT_DIR)/app
EXE_DIR := $(ROOT_DIR)/run
CONFIG_DIR := $(SRC_DIR)/configs
LIBS_DIR := $(ROOT_DIR)/dependencies
LIBS :=
INCS :=

DEFINES :=

EXE := $(EXE_DIR)/$(TARGET)

# Configuration settings
FC := mpifort
FFLAGS :=
AR := ar rcs
LD := $(FC)
RM := rm -f
GD := $(SRC_DIR)/.gen-deps.awk

# edit config/build.conf file desired
include $(ROOT_DIR)/build.conf
include $(CONFIG_DIR)/flags.mk
include $(CONFIG_DIR)/compilers.mk

LIBS += -L$(LIBS_DIR)/2decomp_fft/lib -l2decomp_fft -lfftw3
INCS += -I$(LIBS_DIR)/2decomp_fft/include

# List of all source files
SRCS := $(wildcard $(SRC_DIR)/*.f90) $(wildcard $(APP_DIR)/*.f90)
TEST_SRCS := 

# Define a map from each file name to its object file
obj = $(src).o
$(foreach src, $(SRCS) $(TEST_SRCS), $(eval $(src) := $(obj)))

# Create lists of the build artefacts in this project
OBJS := $(addsuffix .o, $(SRCS))
DEPS := $(addsuffix .d, $(SRCS))
TEST_OBJS := $(addsuffix .o, $(TEST_SRCS))
TEST_DEPS := $(addsuffix .d, $(TEST_SRCS))
LIB := $(patsubst %, lib%.a, $(NAME))
TEST_EXE := $(patsubst %.f90, %.exe, $(TEST_SRCS))

# Declare all public targets
.PHONY: all clean
#all: $(LIB) $(TEST_EXE) $(EXE)
all: $(TEST_EXE) $(EXE)

# Create the static library from the object files
#$(LIB): #$(OBJS)
#	$(AR) $@ $^

# Link the test executables
$(TEST_EXE): %.exe: %.f90.o $(LIB)
	$(LD) -o $@ $^

$(EXE): $(OBJS)
	@mkdir -p $(EXE_DIR)/data
	@cp $(SRC_DIR)/dns.in $(EXE_DIR)
	$(FC) $(FFLAGS) $^ $(LIBS) $(INCS) -o $(EXE)

# Create object files from Fortran source
$(OBJS) $(TEST_OBJS): %.o: % | %.d
	$(FC) $(FFLAGS) -cpp $(DEFINES) $(INCS) $(FFLAGS_MOD_DIR) $(SRC_DIR) -c -o $@ $<

# Process the Fortran source for module dependencies
$(DEPS) $(TEST_DEPS): %.d: %
	$(GD) $< > $@

# Define all module interdependencies
include $(DEPS) $(TEST_DEPS)
$(foreach dep, $(OBJS) $(TEST_OBJS), $(eval $(dep): $($(dep))))

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS) $(TEST_OBJS)) $(filter %.d, $(DEPS) $(TEST_DEPS)) $(filter %.exe, $(TEST_EXE)) $(wildcard *.mod) $(EXE)
