###############################################################################
# FILE: Makefile.am
# AUTHOR: Paul Pavlidis
# CREATE DATE: ~4/2004
# PROJECT: GIST
# DESCRIPTION: Template for Makefile.
###############################################################################

# The programs that comprise Gist.
bin_PROGRAMS = gist-train-svm \
	gist-fselect \
	gist-fast-classify \
	gist-classify \
	gist-matrix \
	gist-project \
	gist-kpca \
	gist-kfd \
	gist-kernel

# Define a common library for use by Gist.
noinst_LTLIBRARIES = libcommon.la

# Put all of these source files into the common library.
libcommon_la_SOURCES = 	\
	array.h \
	array.c \
	class-array.h \
	class-array.c \
	classify.h \
	classify.c \
	cmdline.h \
	compute-kernel.h \
	compute-kernel.c \
	compute-weights.h \
	compute-weights.c \
	eigens.h \
	eigens.c \
	fast-classify.c \
	kernel.h \
	kernel.c \
	linear-algebra.h \
	linear-algebra.c \
	matrix.h \
	matrix.c \
	project.h \
	project.c \
	rdb-matrix.h \
	rdb-matrix.c \
	string-list.h \
	string-list.c \
	test-main.c \
	test-main.h \
	utils.h \
	utils.c

# Not sure what this does.
gist_train_svm_SOURCES = train-main.c fselect.c fselect.h rdb-matrix.c rdb-matrix.h 
gist_fselect_SOURCES = fselect.c fselect.h rdb-matrix.c rdb-matrix.h 
gist_classify_SOURCES = test-main.c rdb-matrix.c rdb-matrix.h fselect.c fselect.h
gist_fast_classify_SOURCES = fast-classify.c test-main.c rdb-matrix.c rdb-matrix.h fselect.c fselect.h
gist_project_SOURCES = test-main.c rdb-matrix.c rdb-matrix.h project.c project.h  fselect.c fselect.h
gist_kpca_SOURCES = train-main.c rdb-matrix.c rdb-matrix.h kernel-pca.h kernel-pca.c fselect.c fselect.h
gist_matrix_SOURCES = rdb-matrix.c rdb-matrix.h 
gist_kfd_SOURCES = train-main.c kfd.c kfd.h rdb-matrix.c rdb-matrix.h fselect.c fselect.h
gist_kernel_SOURCES = gist-kernel.c

# Every compiled program links to the common library.
gist_train_svm_LDADD = libcommon.la
gist_fselect_LDADD = libcommon.la
gist_fast_classify_LDADD = libcommon.la
gist_classify_LDADD = libcommon.la
gist_matrix_LDADD = libcommon.la
gist_project_LDADD = libcommon.la
gist_kpca_LDADD = libcommon.la
gist_kfd_LDADD = libcommon.la
gist_kernel_LDADD = libcommon.la

# Use different compiler flags when debugging.
if DEBUG
AM_CFLAGS = -g -O0 -Wall -W -Werror -DDEBUG -DBOUNDS_CHECK -DMYRAND -DVERIFY_EIGENS
else
AM_CFLAGS = -O3 -Wall -W -Werror -DNDEBUG
endif

# Specify compiler flags for linking programs.
gist_train_svm_CFLAGS =  -DTRAIN_MAIN -DCOMPUTE_WEIGHTS $(AM_CFLAGS)
gist_fselect_CFLAGS =  -DFSELECT_MAIN $(AM_CFLAGS)
gist_fast_classify_CFLAGS = $(AM_CFLAGS)
gist_classify_CFLAGS =   -DTEST_MAIN -DCLASSIFY $(AM_CFLAGS)
gist_matrix_CFLAGS =  -DRDB_MATRIX_MAIN $(AM_CFLAGS)
gist_kpca_CFLAGS =   -DTRAIN_MAIN -DKERNEL_PCA $(AM_CFLAGS)
gist_project_CFLAGS =   -DTEST_MAIN -DPROJECT $(AM_CFLAGS)
gist_kfd_CFLAGS =   -DTRAIN_MAIN -DKFD $(AM_CFLAGS)
gist_kernel_CFLAGS = $(AM_CFLAGS)

# Arguments for the 'etags' program.
ETAGS_ARGS = ./*.h ./*.c 
