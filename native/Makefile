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
all: alula-panel alula-euler libalula.a

SUITESPARSE := ../third_party/klu/source
KLU_PARTS := KLU AMD BTF COLAMD
KLU_CPPFLAGS := $(addprefix -I$(SUITESPARSE)/,$(addsuffix /Include,$(KLU_PARTS))) -I$(SUITESPARSE)/SuiteSparse_config
KLU_SRC := $(filter-out %/amd_global.c,$(foreach part,$(KLU_PARTS),$(wildcard $(SUITESPARSE)/$(part)/Source/*.c))) \
           $(SUITESPARSE)/SuiteSparse_config/SuiteSparse_config.c
KLU_OBJ := $(patsubst $(SUITESPARSE)/%.c,.klu/%.o,$(KLU_SRC))

.klu/%.o: $(SUITESPARSE)/%.c
	@mkdir -p $(@D)
	$(CC) $(KLU_CPPFLAGS) $(filter-out -Wextra -Wpedantic,$(CFLAGS)) -DNDEBUG -DNPRINT= -DNTIMER -DDINT -c $< -o $@

libalula.a: panel.o panel_streamfunction.o panel_field.o panel_guide.o contour.o trace.o streamtube.o sparse.o channel.o body.o body_jacobian.o body_newton.o body_chart.o airfoil_grid.o grid_relax.o bl.o bl_march.o airfoil_bl.o airfoil_coupled.o airfoil_loads.o $(KLU_OBJ)
	$(AR) rcs $@ $^

panel.o panel_streamfunction.o panel_field.o check-panel-field.o: panel.h panel_internal.h
cli.o cli_euler.o cli_geometry.o: cli_geometry.h panel.h
cli_euler.o: airfoil_grid.h airfoil_bl.h airfoil_coupled.h bl.h body_newton.h body.h
trace.o check-panel-field.o: trace.h panel.h
body_chart.o contour.o: contour.h panel.h
contour.o: panel_internal.h
panel_guide.o: panel_guide.h panel_internal.h contour.h trace.h panel.h
airfoil_grid.o: airfoil_grid_internal.h airfoil_grid.h body_chart.h body_newton.h body.h panel_guide.h panel_internal.h contour.h trace.h panel.h streamtube_internal.h streamtube.h
airfoil_bl.o: airfoil_bl.h airfoil_grid_internal.h airfoil_grid.h body_chart.h body_newton.h body.h body_internal.h contour.h panel.h streamtube_internal.h streamtube.h bl.h
airfoil_coupled.o: airfoil_coupled.h airfoil_bl.h airfoil_grid_internal.h airfoil_grid.h body_chart.h body_newton.h body.h body_internal.h contour.h panel.h bl.h
airfoil_loads.o: airfoil_grid_internal.h airfoil_grid.h body_chart.h body_newton.h body.h body_internal.h contour.h panel.h
airfoil_grid.o grid_relax.o: grid_relax.h
grid_relax.o: panel_internal.h panel.h
bl.o bl_march.o check-bl.o: bl.h
check-panel-field.o: panel_guide.h contour.h
streamtube.o channel.o body.o body_jacobian.o: streamtube_internal.h
streamtube.o channel.o body.o body_jacobian.o check-streamtube.o: streamtube.h panel.h
sparse.o: CPPFLAGS += $(KLU_CPPFLAGS)
sparse.o channel.o: sparse.h
channel.o check-channel.o: channel.h sparse.h
body.o body_jacobian.o body_newton.o body_chart.o: body_internal.h body.h streamtube.h streamtube_internal.h
body_newton.o body_chart.o check-body.o: body_newton.h sparse.h
body_chart.o check-body.o: body_chart.h
body.o body_jacobian.o body_newton.o check-body.o: body.h

alula-panel: cli.o cli_geometry.o libalula.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

alula-euler: cli_euler.o cli_geometry.o libalula.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

check-streamtube: check-streamtube.o libalula.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

check-channel: check-channel.o libalula.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

check-body: check-body.o libalula.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

check-panel-field: check-panel-field.o libalula.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

check-bl: check-bl.o libalula.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

check: alula-panel alula-euler check-panel-field check-streamtube check-channel check-body check-bl
	node check-panel.mjs
	node check-streamtube.mjs
	node check-channel.mjs
	node check-body.mjs
	node check-airfoil.mjs
	node check-bl.mjs

clean:
	rm -f alula-panel alula-euler check-panel-field check-streamtube check-channel check-body check-bl libalula.a *.o
	rm -rf .klu
