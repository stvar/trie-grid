# Copyright (C) 2016  Stefan Vargyas
# 
# This file is part of Trie-Grid.
# 
# Trie-Grid is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Trie-Grid is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Trie-Grid.  If not, see <http://www.gnu.org/licenses/>.

PROG := trie

ifdef PROF
PROG := ${PROG}-profiled
endif

ifdef OPT
PROG := ${PROG}-optimized
endif

GCC := gcc
GCC_STD := gnu++98

GCC_VERSION := $(shell \
    bash -c '. ./trie.sh && trie-gcc-version ${GCC}' 2>/dev/null)

ifeq (${GCC_VERSION},)
$(error cannot determine the version of ${GCC})
endif

GCC_SUPPORTS_CXX11 := $(shell \
    bash -c 'test "${GCC_VERSION}" -ge 40800 && echo 1 || echo 0' 2>/dev/null)

GCC_GENERATES_DWARF4 := $(shell \
    bash -c 'test "${GCC_VERSION}" -ge 40800 && echo 1 || echo 0' 2>/dev/null)

ifeq (${GCC_SUPPORTS_CXX11},)
$(error cannot determine whether ${GCC} of version ${GCC_VERSION} supports C++11)
endif

ifeq (${GCC_SUPPORTS_CXX11},0)
ifeq (${GCC_STD},gnu++11)
$(error ${GCC} of version ${GCC_VERSION} doesn't support '-std=${GCC_STD}')
endif
else \
ifeq (${GCC_SUPPORTS_CXX11},1)
ifeq (${GCC_STD},)
GCC_STD := gnu++11
endif
endif

CFLAGS := -std=${GCC_STD} \
          -Wall -Wextra -Wno-ignored-qualifiers

ifndef OPT
CFLAGS += -DDEBUG
ifdef DEBUG
CFLAGS += $(addprefix -DDEBUG_,${DEBUG})
endif
endif

ifdef CONFIG
CFLAGS += $(addprefix -DCONFIG_,${CONFIG})
endif

CFLAGS += -DPROGRAM=${PROG}

ifdef PROF
CFLAGS += -pg
ifdef OPT
CFLAGS += -g
ifeq (${GCC_GENERATES_DWARF4},1)
CFLAGS += -gdwarf-3
endif
endif
endif

ifdef OPT
CFLAGS += -O${OPT}
else
CFLAGS += -g
ifeq (${GCC_GENERATES_DWARF4},1)
CFLAGS += -gdwarf-3
endif
endif

LIBS := -lstdc++

default: ${PROG}

${PROG}: trie.cpp
	${GCC} ${CFLAGS} $< -o $@ ${LIBS}


