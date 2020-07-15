TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.c \
    main.c \
    ../model.cpp \
    ../SimulationTemplate.cpp \
    ../ssa.cpp \
    ../VariableSimulationTemplate.cpp \
    ../tau.cpp \
    ../tausimulationtemplate.cpp

SUBDIRS += \
    TauLeapCSolver.pro

HEADERS += \
    ../model.h \
    ../ssa.h \
    ../tau.h

DISTFILES += \
    ../makefile
