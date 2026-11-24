CC ?= cc
CFLAGS ?= -O3
CFLAGS += -std=c11 -Wall -Wextra -Wpedantic
AR ?= ar

# Select either a dependency-free build or any LP64 LAPACK implementation.
LAPACK_LIBS ?=
ifneq ($(strip $(LAPACK_LIBS)),)
CPPFLAGS += -DALULA_LAPACK
LDLIBS += $(LAPACK_LIBS)
endif
ifeq ($(OPENMP),1)
CFLAGS += -fopenmp
LDFLAGS += -fopenmp
endif
LDLIBS += -lm

.PHONY: all clean check
all: alula-panel libalula.a

SUITESPARSE := ../third_party/klu/source
KLU_PARTS := KLU AMD BTF COLAMD
KLU_CPPFLAGS := $(addprefix -I$(SUITESPARSE)/,$(addsuffix /Include,$(KLU_PARTS))) -I$(SUITESPARSE)/SuiteSparse_config
KLU_SRC := $(filter-out %/amd_global.c,$(foreach part,$(KLU_PARTS),$(wildcard $(SUITESPARSE)/$(part)/Source/*.c))) \
           $(SUITESPARSE)/SuiteSparse_config/SuiteSparse_config.c
KLU_OBJ := $(patsubst $(SUITESPARSE)/%.c,.klu/%.o,$(KLU_SRC))

.klu/%.o: $(SUITESPARSE)/%.c
	@mkdir -p $(@D)
	$(CC) $(KLU_CPPFLAGS) $(filter-out -Wextra -Wpedantic,$(CFLAGS)) -DNDEBUG -DNPRINT= -DNTIMER -DDINT -c $< -o $@

libalula.a: panel.o panel_streamfunction.o streamtube.o sparse.o channel.o body.o body_jacobian.o body_newton.o body_chart.o $(KLU_OBJ)
	$(AR) rcs $@ $^

panel.o panel_streamfunction.o: panel.h panel_internal.h
cli.o: panel.h
streamtube.o channel.o body.o body_jacobian.o: streamtube_internal.h
streamtube.o channel.o body.o body_jacobian.o check-streamtube.o: streamtube.h panel.h
sparse.o: CPPFLAGS += $(KLU_CPPFLAGS)
sparse.o channel.o: sparse.h
channel.o check-channel.o: channel.h sparse.h
body.o body_jacobian.o body_newton.o body_chart.o: body_internal.h body.h streamtube.h streamtube_internal.h
body_newton.o body_chart.o check-body.o: body_newton.h sparse.h
body_chart.o check-body.o: body_chart.h
body.o body_jacobian.o body_newton.o check-body.o: body.h

alula-panel: cli.o libalula.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

check-streamtube: check-streamtube.o libalula.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

check-channel: check-channel.o libalula.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

check-body: check-body.o libalula.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

check: alula-panel check-streamtube check-channel check-body
	node check-panel.mjs
	node check-streamtube.mjs
	node check-channel.mjs
	node check-body.mjs

clean:
	rm -f alula-panel check-streamtube check-channel check-body libalula.a *.o
	rm -rf .klu
