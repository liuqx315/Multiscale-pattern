#!/usr/bin/env python
#------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
#------------------------------------------------------------
# Copyright (c) 2014, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
#------------------------------------------------------------
# script to create the set of 'standard solutions' for 
# regression testing

# imports
import sys
from time import time
import numpy as np
import arkode_tools as ark
import os
import subprocess
import shlex
import pickle

# initialize list of all tests
AllTests = []


# #### debugging tests ####
# p = ark.SolParams(adapt_method=0)
# tests = ( ark.Test('Adapt1:nonlin', './ark_analytic_nonlin.exe', p),
#           ark.Test('Adapt1:analytic', './ark_analytic.exe', p))
# AllTests.append(tests)


#### baseline tests ####
p = ark.SolParams()
tests = ( ark.Test('baseline:analytic', './ark_analytic.exe', p),
          ark.Test('baseline:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('baseline:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('baseline:system', './ark_analytic_sys.exe', p),
          ark.Test('baseline:brusselator', './ark_brusselator.exe', p),
          ark.Test('baseline:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('baseline:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('baseline:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('baseline:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('baseline:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('baseline:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('baseline:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('baseline:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('baseline:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('baseline:bruss', './ark_bruss.exe', p),
          ark.Test('baseline:heat1D', './ark_heat1D.exe', p),
          ark.Test('baseline:heat1D_adapt', './ark_heat1D_adapt.exe', p),
          ark.Test('baseline:heat2D', './ark_heat2D.exe', p),
          ark.Test('baseline:hires', './ark_hires.exe', p),
          ark.Test('baseline:medakzo', './ark_medakzo.exe', p),
          ark.Test('baseline:orego', './ark_orego.exe', p),
          ark.Test('baseline:pollu', './ark_pollu.exe', p),
          ark.Test('baseline:ringmod', './ark_ringmod.exe', p),
          ark.Test('baseline:rober', './ark_rober.exe', p),
          ark.Test('baseline:robertson', './ark_robertson.exe', p),
          ark.Test('baseline:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('baseline:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('baseline:vdpol', './ark_vdpol.exe', p),
          ark.Test('baseline:vdpolm', './ark_vdpolm.exe', p) )
AllTests.append(tests)


#### DIRK table 11 tests ####
p = ark.SolParams(imex=0, btable=11)
p2 = ark.SolParams(imex=0, btable=11, rtol=1.e-6)
tests = ( ark.Test('DIRK_t11:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t11:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t11:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t11:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t11:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t11:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t11:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t11:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t11:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t11:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('DIRK_t11:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('DIRK_t11:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('DIRK_t11:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('DIRK_t11:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          # ark.Test('DIRK_t11:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t11:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t11:heat2D', './ark_heat2D.exe', p),
          # ark.Test('DIRK_t11:hires', './ark_hires.exe', p),
          # ark.Test('DIRK_t11:medakzo', './ark_medakzo.exe', p),
          # ark.Test('DIRK_t11:orego', './ark_orego.exe', p),
          ark.Test('DIRK_t11:pollu', './ark_pollu.exe', p),
          ark.Test('DIRK_t11:ringmod', './ark_ringmod.exe', p),
          # ark.Test('DIRK_t11:rober', './ark_rober.exe', p),
          ark.Test('DIRK_t11:robertson', './ark_robertson.exe', p),
          ark.Test('DIRK_t11:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t11:robertson_root', './ark_robertson_root.exe', p),
          # ark.Test('DIRK_t11:vdpol', './ark_vdpol.exe', p),
          # ark.Test('DIRK_t11:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t11:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t11:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t11:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t11:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t11:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t11:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t11:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t11:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t11:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t11:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('DIRK_t11:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('DIRK_t11:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('DIRK_t11:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('DIRK_t11:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          # ark.Test('DIRK_t11:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t11:tight-heat1D', './ark_heat1D.exe', p2),
          # ark.Test('DIRK_t11:tight-heat2D', './ark_heat2D.exe', p2),
          # ark.Test('DIRK_t11:tight-hires', './ark_hires.exe', p2),
          # ark.Test('DIRK_t11:tight-medakzo', './ark_medakzo.exe', p2),
          # ark.Test('DIRK_t11:tight-orego', './ark_orego.exe', p2),
          ark.Test('DIRK_t11:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('DIRK_t11:tight-ringmod', './ark_ringmod.exe', p2),
          # ark.Test('DIRK_t11:tight-rober', './ark_rober.exe', p2),
          ark.Test('DIRK_t11:tight-robertson', './ark_robertson.exe', p2),
          ark.Test('DIRK_t11:tight-robertson_klu', './ark_robertson_klu.exe', p2) )
          #ark.Test('DIRK_t11:tight-robertson_root', './ark_robertson_root.exe', p2),
          # ark.Test('DIRK_t11:tight-vdpol', './ark_vdpol.exe', p2),
          # ark.Test('DIRK_t11:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### DIRK table 12 tests ####
p = ark.SolParams(imex=0, btable=12)
p2 = ark.SolParams(imex=0, btable=12, rtol=1.e-6)
tests = ( ark.Test('DIRK_t12:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t12:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t12:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t12:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t12:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t12:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t12:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t12:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t12:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t12:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('DIRK_t12:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('DIRK_t12:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          # ark.Test('DIRK_t12:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          # ark.Test('DIRK_t12:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          # ark.Test('DIRK_t12:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t12:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t12:heat2D', './ark_heat2D.exe', p),
          # ark.Test('DIRK_t12:hires', './ark_hires.exe', p),
          ark.Test('DIRK_t12:medakzo', './ark_medakzo.exe', p),
          # ark.Test('DIRK_t12:orego', './ark_orego.exe', p),
          # ark.Test('DIRK_t12:pollu', './ark_pollu.exe', p),
          # ark.Test('DIRK_t12:ringmod', './ark_ringmod.exe', p),
          # ark.Test('DIRK_t12:rober', './ark_rober.exe', p),
          # ark.Test('DIRK_t12:robertson', './ark_robertson.exe', p),
          # ark.Test('DIRK_t12:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t12:robertson_root', './ark_robertson_root.exe', p),
          # ark.Test('DIRK_t12:vdpol', './ark_vdpol.exe', p),
          # ark.Test('DIRK_t12:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t12:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t12:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t12:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t12:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t12:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t12:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t12:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t12:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t12:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t12:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('DIRK_t12:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('DIRK_t12:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          # ark.Test('DIRK_t12:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          # ark.Test('DIRK_t12:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          # ark.Test('DIRK_t12:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t12:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('DIRK_t12:tight-heat2D', './ark_heat2D.exe', p2),
          # ark.Test('DIRK_t12:tight-hires', './ark_hires.exe', p2),
          ark.Test('DIRK_t12:tight-medakzo', './ark_medakzo.exe', p2) )
          # ark.Test('DIRK_t12:tight-orego', './ark_orego.exe', p2),
          # ark.Test('DIRK_t12:tight-pollu', './ark_pollu.exe', p2),
          # ark.Test('DIRK_t12:tight-ringmod', './ark_ringmod.exe', p2),
          # ark.Test('DIRK_t12:tight-rober', './ark_rober.exe', p2),
          # ark.Test('DIRK_t12:tight-robertson', './ark_robertson.exe', p2),
          # ark.Test('DIRK_t12:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('DIRK_t12:tight-robertson_root', './ark_robertson_root.exe', p2),
          # ark.Test('DIRK_t12:tight-vdpol', './ark_vdpol.exe', p2),
          # ark.Test('DIRK_t12:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### DIRK table 13 tests ####
p = ark.SolParams(imex=0, btable=13)
p2 = ark.SolParams(imex=0, btable=13, rtol=1.e-6)
tests = ( ark.Test('DIRK_t13:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t13:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t13:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t13:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t13:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t13:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t13:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t13:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t13:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t13:bruss_Me', './ark_brusselator_Me.exe', p),
          # ark.Test('DIRK_t13:bruss1D', './ark_brusselator1D.exe', p),
          # ark.Test('DIRK_t13:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          # ark.Test('DIRK_t13:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          # ark.Test('DIRK_t13:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          # ark.Test('DIRK_t13:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t13:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t13:heat2D', './ark_heat2D.exe', p),
          # ark.Test('DIRK_t13:hires', './ark_hires.exe', p),
          ark.Test('DIRK_t13:medakzo', './ark_medakzo.exe', p),
          # ark.Test('DIRK_t13:orego', './ark_orego.exe', p),
          # ark.Test('DIRK_t13:pollu', './ark_pollu.exe', p),
          # ark.Test('DIRK_t13:ringmod', './ark_ringmod.exe', p),
          # ark.Test('DIRK_t13:rober', './ark_rober.exe', p),
          # ark.Test('DIRK_t13:robertson', './ark_robertson.exe', p),
          # ark.Test('DIRK_t13:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t13:robertson_root', './ark_robertson_root.exe', p),
          # ark.Test('DIRK_t13:vdpol', './ark_vdpol.exe', p),
          # ark.Test('DIRK_t13:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t13:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t13:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t13:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t13:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t13:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t13:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t13:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t13:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t13:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t13:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          # ark.Test('DIRK_t13:tight-bruss1D', './ark_brusselator1D.exe', p2),
          # ark.Test('DIRK_t13:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          # ark.Test('DIRK_t13:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          # ark.Test('DIRK_t13:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          # ark.Test('DIRK_t13:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t13:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('DIRK_t13:tight-heat2D', './ark_heat2D.exe', p2),
          # ark.Test('DIRK_t13:tight-hires', './ark_hires.exe', p2),
          ark.Test('DIRK_t13:tight-medakzo', './ark_medakzo.exe', p2) )
          # ark.Test('DIRK_t13:tight-orego', './ark_orego.exe', p2),
          # ark.Test('DIRK_t13:tight-pollu', './ark_pollu.exe', p2),
          # ark.Test('DIRK_t13:tight-ringmod', './ark_ringmod.exe', p2),
          # ark.Test('DIRK_t13:tight-rober', './ark_rober.exe', p2),
          # ark.Test('DIRK_t13:tight-robertson', './ark_robertson.exe', p2),
          # ark.Test('DIRK_t13:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('DIRK_t13:tight-robertson_root', './ark_robertson_root.exe', p2),
          # ark.Test('DIRK_t13:tight-vdpol', './ark_vdpol.exe', p2),
          # ark.Test('DIRK_t13:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### DIRK table 14 tests ####
p = ark.SolParams(imex=0, btable=14)
p2 = ark.SolParams(imex=0, btable=14, rtol=1.e-6)
tests = ( ark.Test('DIRK_t14:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t14:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t14:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t14:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t14:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t14:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t14:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t14:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t14:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t14:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('DIRK_t14:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('DIRK_t14:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('DIRK_t14:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('DIRK_t14:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('DIRK_t14:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t14:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t14:heat2D', './ark_heat2D.exe', p),
          ark.Test('DIRK_t14:hires', './ark_hires.exe', p),
          ark.Test('DIRK_t14:medakzo', './ark_medakzo.exe', p),
          ark.Test('DIRK_t14:orego', './ark_orego.exe', p),
          ark.Test('DIRK_t14:pollu', './ark_pollu.exe', p),
          ark.Test('DIRK_t14:ringmod', './ark_ringmod.exe', p),
          ark.Test('DIRK_t14:rober', './ark_rober.exe', p),
          # ark.Test('DIRK_t14:robertson', './ark_robertson.exe', p),
          # ark.Test('DIRK_t14:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t14:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('DIRK_t14:vdpol', './ark_vdpol.exe', p),
          # ark.Test('DIRK_t14:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t14:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t14:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t14:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t14:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t14:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t14:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t14:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t14:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t14:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t14:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('DIRK_t14:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('DIRK_t14:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('DIRK_t14:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('DIRK_t14:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('DIRK_t14:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t14:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('DIRK_t14:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('DIRK_t14:tight-hires', './ark_hires.exe', p2),
          ark.Test('DIRK_t14:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('DIRK_t14:tight-orego', './ark_orego.exe', p2),
          ark.Test('DIRK_t14:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('DIRK_t14:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('DIRK_t14:tight-rober', './ark_rober.exe', p2),
          # ark.Test('DIRK_t14:tight-robertson', './ark_robertson.exe', p2),
          # ark.Test('DIRK_t14:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('DIRK_t14:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('DIRK_t14:tight-vdpol', './ark_vdpol.exe', p2) )
          # ark.Test('DIRK_t14:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### DIRK table 15 tests ####
p = ark.SolParams(imex=0, btable=15)
p2 = ark.SolParams(imex=0, btable=15, rtol=1.e-6)
tests = ( ark.Test('DIRK_t15:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t15:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t15:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t15:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t15:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t15:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t15:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t15:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t15:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t15:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('DIRK_t15:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('DIRK_t15:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('DIRK_t15:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('DIRK_t15:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('DIRK_t15:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t15:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t15:heat2D', './ark_heat2D.exe', p),
          ark.Test('DIRK_t15:hires', './ark_hires.exe', p),
          ark.Test('DIRK_t15:medakzo', './ark_medakzo.exe', p),
          ark.Test('DIRK_t15:orego', './ark_orego.exe', p),
          ark.Test('DIRK_t15:pollu', './ark_pollu.exe', p),
          ark.Test('DIRK_t15:ringmod', './ark_ringmod.exe', p),
          # ark.Test('DIRK_t15:rober', './ark_rober.exe', p),
          # ark.Test('DIRK_t15:robertson', './ark_robertson.exe', p),
          # ark.Test('DIRK_t15:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t15:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('DIRK_t15:vdpol', './ark_vdpol.exe', p),
          # ark.Test('DIRK_t15:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t15:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t15:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t15:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t15:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t15:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t15:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t15:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t15:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t15:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t15:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('DIRK_t15:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('DIRK_t15:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('DIRK_t15:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('DIRK_t15:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('DIRK_t15:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t15:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('DIRK_t15:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('DIRK_t15:tight-hires', './ark_hires.exe', p2),
          ark.Test('DIRK_t15:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('DIRK_t15:tight-orego', './ark_orego.exe', p2),
          ark.Test('DIRK_t15:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('DIRK_t15:tight-ringmod', './ark_ringmod.exe', p2),
          # ark.Test('DIRK_t15:tight-rober', './ark_rober.exe', p2),
          # ark.Test('DIRK_t15:tight-robertson', './ark_robertson.exe', p2),
          # ark.Test('DIRK_t15:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('DIRK_t15:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('DIRK_t15:tight-vdpol', './ark_vdpol.exe', p2) )
          # ark.Test('DIRK_t15:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### DIRK table 16 tests ####
p = ark.SolParams(imex=0, btable=16)
p2 = ark.SolParams(imex=0, btable=16, rtol=1.e-6)
tests = ( ark.Test('DIRK_t16:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t16:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t16:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t16:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t16:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t16:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t16:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t16:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t16:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t16:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('DIRK_t16:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('DIRK_t16:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('DIRK_t16:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('DIRK_t16:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('DIRK_t16:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t16:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t16:heat2D', './ark_heat2D.exe', p),
          # ark.Test('DIRK_t16:hires', './ark_hires.exe', p),
          ark.Test('DIRK_t16:medakzo', './ark_medakzo.exe', p),
          ark.Test('DIRK_t16:orego', './ark_orego.exe', p),
          ark.Test('DIRK_t16:pollu', './ark_pollu.exe', p),
          ark.Test('DIRK_t16:ringmod', './ark_ringmod.exe', p),
          ark.Test('DIRK_t16:rober', './ark_rober.exe', p),
          ark.Test('DIRK_t16:robertson', './ark_robertson.exe', p),
          ark.Test('DIRK_t16:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t16:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('DIRK_t16:vdpol', './ark_vdpol.exe', p),
          ark.Test('DIRK_t16:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t16:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t16:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t16:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t16:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t16:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t16:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t16:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t16:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t16:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t16:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('DIRK_t16:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('DIRK_t16:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('DIRK_t16:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('DIRK_t16:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('DIRK_t16:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t16:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('DIRK_t16:tight-heat2D', './ark_heat2D.exe', p2),
          # ark.Test('DIRK_t16:tight-hires', './ark_hires.exe', p2),
          ark.Test('DIRK_t16:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('DIRK_t16:tight-orego', './ark_orego.exe', p2),
          ark.Test('DIRK_t16:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('DIRK_t16:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('DIRK_t16:tight-rober', './ark_rober.exe', p2),
          ark.Test('DIRK_t16:tight-robertson', './ark_robertson.exe', p2),
          ark.Test('DIRK_t16:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('DIRK_t16:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('DIRK_t16:tight-vdpol', './ark_vdpol.exe', p2),
          ark.Test('DIRK_t16:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### DIRK table 17 tests ####
p = ark.SolParams(imex=0, btable=17)
p2 = ark.SolParams(imex=0, btable=17, rtol=1.e-6)
tests = ( ark.Test('DIRK_t17:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t17:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t17:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t17:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t17:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t17:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t17:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t17:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t17:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t17:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('DIRK_t17:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('DIRK_t17:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('DIRK_t17:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('DIRK_t17:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('DIRK_t17:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t17:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t17:heat2D', './ark_heat2D.exe', p),
          ark.Test('DIRK_t17:hires', './ark_hires.exe', p),
          ark.Test('DIRK_t17:medakzo', './ark_medakzo.exe', p),
          # ark.Test('DIRK_t17:orego', './ark_orego.exe', p),
          ark.Test('DIRK_t17:pollu', './ark_pollu.exe', p),
          ark.Test('DIRK_t17:ringmod', './ark_ringmod.exe', p),
          ark.Test('DIRK_t17:rober', './ark_rober.exe', p),
          ark.Test('DIRK_t17:robertson', './ark_robertson.exe', p),
          ark.Test('DIRK_t17:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t17:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('DIRK_t17:vdpol', './ark_vdpol.exe', p),
          ark.Test('DIRK_t17:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t17:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t17:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t17:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t17:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t17:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t17:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t17:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t17:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t17:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t17:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('DIRK_t17:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('DIRK_t17:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('DIRK_t17:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('DIRK_t17:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('DIRK_t17:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t17:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('DIRK_t17:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('DIRK_t17:tight-hires', './ark_hires.exe', p2),
          ark.Test('DIRK_t17:tight-medakzo', './ark_medakzo.exe', p2),
          # ark.Test('DIRK_t17:tight-orego', './ark_orego.exe', p2),
          ark.Test('DIRK_t17:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('DIRK_t17:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('DIRK_t17:tight-rober', './ark_rober.exe', p2),
          ark.Test('DIRK_t17:tight-robertson', './ark_robertson.exe', p2),
          ark.Test('DIRK_t17:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('DIRK_t17:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('DIRK_t17:tight-vdpol', './ark_vdpol.exe', p2),
          ark.Test('DIRK_t17:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### DIRK table 18 tests ####
p = ark.SolParams(imex=0, btable=18)
p2 = ark.SolParams(imex=0, btable=18, rtol=1.e-6)
tests = ( ark.Test('DIRK_t18:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t18:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t18:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t18:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t18:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t18:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t18:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t18:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t18:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t18:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('DIRK_t18:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('DIRK_t18:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('DIRK_t18:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('DIRK_t18:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('DIRK_t18:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t18:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t18:heat2D', './ark_heat2D.exe', p),
          ark.Test('DIRK_t18:hires', './ark_hires.exe', p),
          ark.Test('DIRK_t18:medakzo', './ark_medakzo.exe', p),
          ark.Test('DIRK_t18:orego', './ark_orego.exe', p),
          ark.Test('DIRK_t18:pollu', './ark_pollu.exe', p),
          ark.Test('DIRK_t18:ringmod', './ark_ringmod.exe', p),
          ark.Test('DIRK_t18:rober', './ark_rober.exe', p),
          ark.Test('DIRK_t18:robertson', './ark_robertson.exe', p),
          ark.Test('DIRK_t18:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t18:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('DIRK_t18:vdpol', './ark_vdpol.exe', p),
          ark.Test('DIRK_t18:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t18:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t18:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t18:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t18:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t18:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t18:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t18:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t18:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t18:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t18:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('DIRK_t18:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('DIRK_t18:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('DIRK_t18:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('DIRK_t18:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('DIRK_t18:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t18:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('DIRK_t18:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('DIRK_t18:tight-hires', './ark_hires.exe', p2),
          ark.Test('DIRK_t18:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('DIRK_t18:tight-orego', './ark_orego.exe', p2),
          ark.Test('DIRK_t18:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('DIRK_t18:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('DIRK_t18:tight-rober', './ark_rober.exe', p2),
          ark.Test('DIRK_t18:tight-robertson', './ark_robertson.exe', p2),
          ark.Test('DIRK_t18:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('DIRK_t18:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('DIRK_t18:tight-vdpol', './ark_vdpol.exe', p2),
          ark.Test('DIRK_t18:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### DIRK table 19 tests ####
p = ark.SolParams(imex=0, btable=19)
p2 = ark.SolParams(imex=0, btable=19, rtol=1.e-6)
tests = ( ark.Test('DIRK_t19:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t19:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t19:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t19:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t19:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t19:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t19:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t19:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t19:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t19:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('DIRK_t19:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('DIRK_t19:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('DIRK_t19:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('DIRK_t19:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('DIRK_t19:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t19:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t19:heat2D', './ark_heat2D.exe', p),
          ark.Test('DIRK_t19:hires', './ark_hires.exe', p),
          ark.Test('DIRK_t19:medakzo', './ark_medakzo.exe', p),
          ark.Test('DIRK_t19:orego', './ark_orego.exe', p),
          ark.Test('DIRK_t19:pollu', './ark_pollu.exe', p),
          ark.Test('DIRK_t19:ringmod', './ark_ringmod.exe', p),
          ark.Test('DIRK_t19:rober', './ark_rober.exe', p),
          # ark.Test('DIRK_t19:robertson', './ark_robertson.exe', p),
          # ark.Test('DIRK_t19:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t19:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('DIRK_t19:vdpol', './ark_vdpol.exe', p),
          # ark.Test('DIRK_t19:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t19:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t19:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t19:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t19:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t19:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t19:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t19:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t19:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t19:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t19:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('DIRK_t19:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('DIRK_t19:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('DIRK_t19:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('DIRK_t19:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('DIRK_t19:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t19:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('DIRK_t19:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('DIRK_t19:tight-hires', './ark_hires.exe', p2),
          ark.Test('DIRK_t19:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('DIRK_t19:tight-orego', './ark_orego.exe', p2),
          ark.Test('DIRK_t19:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('DIRK_t19:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('DIRK_t19:tight-rober', './ark_rober.exe', p2),
          # ark.Test('DIRK_t19:tight-robertson', './ark_robertson.exe', p2),
          # ark.Test('DIRK_t19:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('DIRK_t19:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('DIRK_t19:tight-vdpol', './ark_vdpol.exe', p2) )
          # ark.Test('DIRK_t19:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### DIRK table 20 tests ####
p = ark.SolParams(imex=0, btable=20)
p2 = ark.SolParams(imex=0, btable=20, rtol=1.e-6)
tests = ( ark.Test('DIRK_t20:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t20:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t20:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t20:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t20:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t20:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t20:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t20:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t20:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t20:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('DIRK_t20:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('DIRK_t20:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('DIRK_t20:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('DIRK_t20:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('DIRK_t20:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t20:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t20:heat2D', './ark_heat2D.exe', p),
          ark.Test('DIRK_t20:hires', './ark_hires.exe', p),
          ark.Test('DIRK_t20:medakzo', './ark_medakzo.exe', p),
          # ark.Test('DIRK_t20:orego', './ark_orego.exe', p),
          ark.Test('DIRK_t20:pollu', './ark_pollu.exe', p),
          ark.Test('DIRK_t20:ringmod', './ark_ringmod.exe', p),
          # ark.Test('DIRK_t20:rober', './ark_rober.exe', p),
          # ark.Test('DIRK_t20:robertson', './ark_robertson.exe', p),
          # ark.Test('DIRK_t20:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t20:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('DIRK_t20:vdpol', './ark_vdpol.exe', p),
          # ark.Test('DIRK_t20:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t20:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t20:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t20:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t20:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t20:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t20:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t20:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t20:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t20:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t20:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('DIRK_t20:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('DIRK_t20:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('DIRK_t20:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('DIRK_t20:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('DIRK_t20:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t20:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('DIRK_t20:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('DIRK_t20:tight-hires', './ark_hires.exe', p2),
          ark.Test('DIRK_t20:tight-medakzo', './ark_medakzo.exe', p2),
          # ark.Test('DIRK_t20:tight-orego', './ark_orego.exe', p2),
          ark.Test('DIRK_t20:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('DIRK_t20:tight-ringmod', './ark_ringmod.exe', p2),
          # ark.Test('DIRK_t20:tight-rober', './ark_rober.exe', p2),
          # ark.Test('DIRK_t20:tight-robertson', './ark_robertson.exe', p2),
          # ark.Test('DIRK_t20:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('DIRK_t20:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('DIRK_t20:tight-vdpol', './ark_vdpol.exe', p2) )
          # ark.Test('DIRK_t20:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### DIRK table 21 tests ####
p = ark.SolParams(imex=0, btable=21)
p2 = ark.SolParams(imex=0, btable=21, rtol=1.e-6)
tests = ( ark.Test('DIRK_t21:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t21:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t21:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t21:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t21:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t21:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t21:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t21:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t21:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t21:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('DIRK_t21:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('DIRK_t21:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('DIRK_t21:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('DIRK_t21:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('DIRK_t21:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t21:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t21:heat2D', './ark_heat2D.exe', p),
          ark.Test('DIRK_t21:hires', './ark_hires.exe', p),
          ark.Test('DIRK_t21:medakzo', './ark_medakzo.exe', p),
          ark.Test('DIRK_t21:orego', './ark_orego.exe', p),
          ark.Test('DIRK_t21:pollu', './ark_pollu.exe', p),
          ark.Test('DIRK_t21:ringmod', './ark_ringmod.exe', p),
          # ark.Test('DIRK_t21:rober', './ark_rober.exe', p),
          # ark.Test('DIRK_t21:robertson', './ark_robertson.exe', p),
          # ark.Test('DIRK_t21:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t21:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('DIRK_t21:vdpol', './ark_vdpol.exe', p),
          # ark.Test('DIRK_t21:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t21:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t21:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t21:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t21:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t21:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t21:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t21:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t21:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t21:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t21:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('DIRK_t21:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('DIRK_t21:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('DIRK_t21:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('DIRK_t21:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('DIRK_t21:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t21:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('DIRK_t21:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('DIRK_t21:tight-hires', './ark_hires.exe', p2),
          ark.Test('DIRK_t21:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('DIRK_t21:tight-orego', './ark_orego.exe', p2),
          ark.Test('DIRK_t21:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('DIRK_t21:tight-ringmod', './ark_ringmod.exe', p2),
          # ark.Test('DIRK_t21:tight-rober', './ark_rober.exe', p2),
          # ark.Test('DIRK_t21:tight-robertson', './ark_robertson.exe', p2),
          # ark.Test('DIRK_t21:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('DIRK_t21:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('DIRK_t21:tight-vdpol', './ark_vdpol.exe', p2) )
          # ark.Test('DIRK_t21:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### DIRK table 22 tests ####
p = ark.SolParams(imex=0, btable=22)
p2 = ark.SolParams(imex=0, btable=22, rtol=1.e-6)
tests = ( ark.Test('DIRK_t22:analytic', './ark_analytic.exe', p),
          ark.Test('DIRK_t22:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('DIRK_t22:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('DIRK_t22:system', './ark_analytic_sys.exe', p),
          ark.Test('DIRK_t22:brusselator', './ark_brusselator.exe', p),
          ark.Test('DIRK_t22:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('DIRK_t22:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('DIRK_t22:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('DIRK_t22:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('DIRK_t22:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('DIRK_t22:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('DIRK_t22:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('DIRK_t22:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('DIRK_t22:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('DIRK_t22:bruss', './ark_bruss.exe', p),
          ark.Test('DIRK_t22:heat1D', './ark_heat1D.exe', p),
          ark.Test('DIRK_t22:heat2D', './ark_heat2D.exe', p),
          ark.Test('DIRK_t22:hires', './ark_hires.exe', p),
          ark.Test('DIRK_t22:medakzo', './ark_medakzo.exe', p),
          ark.Test('DIRK_t22:orego', './ark_orego.exe', p),
          ark.Test('DIRK_t22:pollu', './ark_pollu.exe', p),
          ark.Test('DIRK_t22:ringmod', './ark_ringmod.exe', p),
          ark.Test('DIRK_t22:rober', './ark_rober.exe', p),
          ark.Test('DIRK_t22:robertson', './ark_robertson.exe', p),
          ark.Test('DIRK_t22:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('DIRK_t22:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('DIRK_t22:vdpol', './ark_vdpol.exe', p),
          ark.Test('DIRK_t22:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('DIRK_t22:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('DIRK_t22:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('DIRK_t22:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('DIRK_t22:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('DIRK_t22:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('DIRK_t22:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('DIRK_t22:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('DIRK_t22:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('DIRK_t22:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('DIRK_t22:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('DIRK_t22:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('DIRK_t22:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('DIRK_t22:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('DIRK_t22:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('DIRK_t22:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('DIRK_t22:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('DIRK_t22:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('DIRK_t22:tight-hires', './ark_hires.exe', p2),
          ark.Test('DIRK_t22:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('DIRK_t22:tight-orego', './ark_orego.exe', p2),
          ark.Test('DIRK_t22:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('DIRK_t22:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('DIRK_t22:tight-rober', './ark_rober.exe', p2),
          ark.Test('DIRK_t22:tight-robertson', './ark_robertson.exe', p2),
          ark.Test('DIRK_t22:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('DIRK_t22:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('DIRK_t22:tight-vdpol', './ark_vdpol.exe', p2),
          ark.Test('DIRK_t22:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### ARK order 3 tests ####
p = ark.SolParams(imex=2, order=3)
p2 = ark.SolParams(imex=2, order=3, rtol=1.e-6)
tests = ( ark.Test('ARK_o3:analytic', './ark_analytic.exe', p),
          ark.Test('ARK_o3:system', './ark_analytic_sys.exe', p),
          ark.Test('ARK_o3:brusselator', './ark_brusselator.exe', p),
          ark.Test('ARK_o3:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('ARK_o3:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('ARK_o3:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('ARK_o3:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('ARK_o3:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('ARK_o3:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('ARK_o3:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          #ark.Test('ARK_o3:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          #ark.Test('ARK_o3:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('ARK_o3:bruss', './ark_bruss.exe', p),
          ark.Test('ARK_o3:medakzo', './ark_medakzo.exe', p),
          ark.Test('ARK_o3:pollu', './ark_pollu.exe', p),
          ark.Test('ARK_o3:vdpol', './ark_vdpol.exe', p),
          ark.Test('ARK_o3:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('ARK_o3:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ARK_o3:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ARK_o3:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ARK_o3:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('ARK_o3:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('ARK_o3:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('ARK_o3:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('ARK_o3:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('ARK_o3:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('ARK_o3:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          #ark.Test('ARK_o3:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          #ark.Test('ARK_o3:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('ARK_o3:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('ARK_o3:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('ARK_o3:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('ARK_o3:tight-vdpol', './ark_vdpol.exe', p2),
          ark.Test('ARK_o3:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### ARK order 4 tests ####
p = ark.SolParams(imex=2, order=4)
p2 = ark.SolParams(imex=2, order=4, rtol=1.e-6)
tests = ( ark.Test('ARK_o4:analytic', './ark_analytic.exe', p),
          ark.Test('ARK_o4:system', './ark_analytic_sys.exe', p),
          ark.Test('ARK_o4:brusselator', './ark_brusselator.exe', p),
          ark.Test('ARK_o4:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('ARK_o4:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('ARK_o4:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('ARK_o4:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('ARK_o4:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('ARK_o4:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('ARK_o4:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          #ark.Test('ARK_o4:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          #ark.Test('ARK_o4:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('ARK_o4:bruss', './ark_bruss.exe', p),
          ark.Test('ARK_o4:medakzo', './ark_medakzo.exe', p),
          ark.Test('ARK_o4:pollu', './ark_pollu.exe', p),
          ark.Test('ARK_o4:vdpol', './ark_vdpol.exe', p),
          ark.Test('ARK_o4:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ARK_o4:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ARK_o4:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ARK_o4:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('ARK_o4:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('ARK_o4:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('ARK_o4:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('ARK_o4:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('ARK_o4:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('ARK_o4:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          #ark.Test('ARK_o4:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          #ark.Test('ARK_o4:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('ARK_o4:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('ARK_o4:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('ARK_o4:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('ARK_o4:tight-vdpol', './ark_vdpol.exe', p2) )
AllTests.append(tests)


#### ARK order 5 tests ####
p = ark.SolParams(imex=2, order=5)
p2 = ark.SolParams(imex=2, order=5, rtol=1.e-6)
tests = ( ark.Test('ARK_o5:analytic', './ark_analytic.exe', p),
          ark.Test('ARK_o5:system', './ark_analytic_sys.exe', p),
          ark.Test('ARK_o5:brusselator', './ark_brusselator.exe', p),
          ark.Test('ARK_o5:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('ARK_o5:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('ARK_o5:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('ARK_o5:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('ARK_o5:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('ARK_o5:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('ARK_o5:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          #ark.Test('ARK_o5:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          #ark.Test('ARK_o5:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('ARK_o5:bruss', './ark_bruss.exe', p),
          ark.Test('ARK_o5:medakzo', './ark_medakzo.exe', p),
          ark.Test('ARK_o5:pollu', './ark_pollu.exe', p),
          ark.Test('ARK_o5:vdpol', './ark_vdpol.exe', p),
          ark.Test('ARK_o5:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('ARK_o5:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ARK_o5:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ARK_o5:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ARK_o5:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('ARK_o5:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('ARK_o5:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('ARK_o5:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('ARK_o5:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('ARK_o5:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('ARK_o5:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          #ark.Test('ARK_o5:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          #ark.Test('ARK_o5:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('ARK_o5:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('ARK_o5:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('ARK_o5:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('ARK_o5:tight-vdpol', './ark_vdpol.exe', p2),
          ark.Test('ARK_o5:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### ERK table 0 tests ####
p = ark.SolParams(imex=1, btable=0)
p2 = ark.SolParams(imex=1, btable=0, rtol=1.e-6)
tests = ( ark.Test('ERK_t0:analytic', './ark_analytic.exe', p),
          ark.Test('ERK_t0:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('ERK_t0:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('ERK_t0:system', './ark_analytic_sys.exe', p),
          ark.Test('ERK_t0:brusselator', './ark_brusselator.exe', p),
          ark.Test('ERK_t0:heat1D', './ark_heat1D.exe' , p),
          ark.Test('ERK_t0:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ERK_t0:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('ERK_t0:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('ERK_t0:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ERK_t0:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ERK_t0:tight-heat1D', './ark_heat1D.exe' , p2) )
AllTests.append(tests)


#### ERK table 1 tests ####
p = ark.SolParams(imex=1, btable=1)
p2 = ark.SolParams(imex=1, btable=1, rtol=1.e-6)
tests = ( ark.Test('ERK_t1:analytic', './ark_analytic.exe', p),
          ark.Test('ERK_t1:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('ERK_t1:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('ERK_t1:system', './ark_analytic_sys.exe', p),
          ark.Test('ERK_t1:brusselator', './ark_brusselator.exe', p),
          ark.Test('ERK_t1:heat1D', './ark_heat1D.exe' , p),
          ark.Test('ERK_t1:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ERK_t1:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('ERK_t1:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('ERK_t1:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ERK_t1:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ERK_t1:tight-heat1D', './ark_heat1D.exe' , p2) )
AllTests.append(tests)


#### ERK table 2 tests ####
p = ark.SolParams(imex=1, btable=2)
p2 = ark.SolParams(imex=1, btable=2, rtol=1.e-6)
tests = ( ark.Test('ERK_t2:analytic', './ark_analytic.exe', p),
          ark.Test('ERK_t2:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('ERK_t2:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('ERK_t2:system', './ark_analytic_sys.exe', p),
          ark.Test('ERK_t2:brusselator', './ark_brusselator.exe', p),
          ark.Test('ERK_t2:heat1D', './ark_heat1D.exe' , p),
          ark.Test('ERK_t2:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ERK_t2:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('ERK_t2:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('ERK_t2:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ERK_t2:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ERK_t2:tight-heat1D', './ark_heat1D.exe' , p2) )
AllTests.append(tests)


#### ERK table 3 tests ####
p = ark.SolParams(imex=1, btable=3)
p2 = ark.SolParams(imex=1, btable=3, rtol=1.e-6)
tests = ( ark.Test('ERK_t3:analytic', './ark_analytic.exe', p),
          ark.Test('ERK_t3:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('ERK_t3:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('ERK_t3:system', './ark_analytic_sys.exe', p),
          ark.Test('ERK_t3:brusselator', './ark_brusselator.exe', p),
          ark.Test('ERK_t3:heat1D', './ark_heat1D.exe' , p),
          ark.Test('ERK_t3:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ERK_t3:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('ERK_t3:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('ERK_t3:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ERK_t3:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ERK_t3:tight-heat1D', './ark_heat1D.exe' , p2) )
AllTests.append(tests)


#### ERK table 4 tests ####
p = ark.SolParams(imex=1, btable=4)
p2 = ark.SolParams(imex=1, btable=4, rtol=1.e-6)
tests = ( ark.Test('ERK_t4:analytic', './ark_analytic.exe', p),
          ark.Test('ERK_t4:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('ERK_t4:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('ERK_t4:system', './ark_analytic_sys.exe', p),
          ark.Test('ERK_t4:brusselator', './ark_brusselator.exe', p),
          ark.Test('ERK_t4:heat1D', './ark_heat1D.exe' , p),
          ark.Test('ERK_t4:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ERK_t4:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('ERK_t4:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('ERK_t4:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ERK_t4:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ERK_t4:tight-heat1D', './ark_heat1D.exe' , p2) )
AllTests.append(tests)


#### ERK table 5 tests ####
p = ark.SolParams(imex=1, btable=5)
p2 = ark.SolParams(imex=1, btable=5, rtol=1.e-6)
tests = ( ark.Test('ERK_t5:analytic', './ark_analytic.exe', p),
          ark.Test('ERK_t5:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('ERK_t5:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('ERK_t5:system', './ark_analytic_sys.exe', p),
          ark.Test('ERK_t5:brusselator', './ark_brusselator.exe', p),
          ark.Test('ERK_t5:heat1D', './ark_heat1D.exe' , p),
          ark.Test('ERK_t5:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ERK_t5:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('ERK_t5:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('ERK_t5:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ERK_t5:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ERK_t5:tight-heat1D', './ark_heat1D.exe' , p2) )
AllTests.append(tests)


#### ERK table 6 tests ####
p = ark.SolParams(imex=1, btable=6)
p2 = ark.SolParams(imex=1, btable=6, rtol=1.e-6)
tests = ( ark.Test('ERK_t6:analytic', './ark_analytic.exe', p),
          ark.Test('ERK_t6:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('ERK_t6:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('ERK_t6:system', './ark_analytic_sys.exe', p),
          ark.Test('ERK_t6:brusselator', './ark_brusselator.exe', p),
          ark.Test('ERK_t6:heat1D', './ark_heat1D.exe' , p),
          ark.Test('ERK_t6:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ERK_t6:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('ERK_t6:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('ERK_t6:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ERK_t6:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ERK_t6:tight-heat1D', './ark_heat1D.exe' , p2) )
AllTests.append(tests)


#### ERK table 7 tests ####
p = ark.SolParams(imex=1, btable=7)
p2 = ark.SolParams(imex=1, btable=7, rtol=1.e-6)
tests = ( ark.Test('ERK_t7:analytic', './ark_analytic.exe', p),
          ark.Test('ERK_t7:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('ERK_t7:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('ERK_t7:system', './ark_analytic_sys.exe', p),
          ark.Test('ERK_t7:brusselator', './ark_brusselator.exe', p),
          ark.Test('ERK_t7:heat1D', './ark_heat1D.exe' , p),
          ark.Test('ERK_t7:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ERK_t7:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('ERK_t7:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('ERK_t7:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ERK_t7:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ERK_t7:tight-heat1D', './ark_heat1D.exe' , p2) )
AllTests.append(tests)


#### ERK table 8 tests ####
p = ark.SolParams(imex=1, btable=8)
p2 = ark.SolParams(imex=1, btable=8, rtol=1.e-6)
tests = ( ark.Test('ERK_t8:analytic', './ark_analytic.exe', p),
          ark.Test('ERK_t8:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('ERK_t8:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('ERK_t8:system', './ark_analytic_sys.exe', p),
          ark.Test('ERK_t8:brusselator', './ark_brusselator.exe', p),
          ark.Test('ERK_t8:heat1D', './ark_heat1D.exe' , p),
          ark.Test('ERK_t8:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ERK_t8:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('ERK_t8:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('ERK_t8:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ERK_t8:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ERK_t8:tight-heat1D', './ark_heat1D.exe' , p2) )
AllTests.append(tests)


#### ERK table 9 tests ####
p = ark.SolParams(imex=1, btable=9)
p2 = ark.SolParams(imex=1, btable=9, rtol=1.e-6)
tests = ( ark.Test('ERK_t9:analytic', './ark_analytic.exe', p),
          ark.Test('ERK_t9:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('ERK_t9:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('ERK_t9:system', './ark_analytic_sys.exe', p),
          ark.Test('ERK_t9:brusselator', './ark_brusselator.exe', p),
          ark.Test('ERK_t9:heat1D', './ark_heat1D.exe' , p),
          ark.Test('ERK_t9:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ERK_t9:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('ERK_t9:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('ERK_t9:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ERK_t9:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ERK_t9:tight-heat1D', './ark_heat1D.exe' , p2) )
AllTests.append(tests)


#### ERK table 10 tests ####
p = ark.SolParams(imex=1, btable=10)
p2 = ark.SolParams(imex=1, btable=10, rtol=1.e-6)
tests = ( ark.Test('ERK_t10:analytic', './ark_analytic.exe', p),
          ark.Test('ERK_t10:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('ERK_t10:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('ERK_t10:system', './ark_analytic_sys.exe', p),
          ark.Test('ERK_t10:brusselator', './ark_brusselator.exe', p),
          ark.Test('ERK_t10:heat1D', './ark_heat1D.exe' , p),
          ark.Test('ERK_t10:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('ERK_t10:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('ERK_t10:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('ERK_t10:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('ERK_t10:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('ERK_t10:tight-heat1D', './ark_heat1D.exe' , p2) )
AllTests.append(tests)


#### fixed-point solver tests ####
p = ark.SolParams(imex=0, order=3, fixedpt=1, m_aa=3, maxcor=50)
tests = ( ark.Test('fixed:point-analytic', './ark_analytic.exe', p),
          ark.Test('fixed:point-nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('fixed:point-nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('fixed:point-system', './ark_analytic_sys.exe', p),
          ark.Test('fixed:point-brusselator', './ark_brusselator.exe', p),
          ark.Test('fixed:point-bruss', './ark_bruss.exe', p),
          ark.Test('fixed:point-heat1D', './ark_heat1D.exe', p),
          ark.Test('fixed:point-hires', './ark_hires.exe', p),
          ark.Test('fixed:point-medakzo', './ark_medakzo.exe', p) )
AllTests.append(tests)


#### adaptivity method 0 tests ####
p = ark.SolParams(adapt_method=0)
p2 = ark.SolParams(adapt_method=0, rtol=1.e-6)
tests = ( ark.Test('Adapt0:analytic', './ark_analytic.exe', p),
          ark.Test('Adapt0:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('Adapt0:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('Adapt0:system', './ark_analytic_sys.exe', p),
          ark.Test('Adapt0:brusselator', './ark_brusselator.exe', p),
          ark.Test('Adapt0:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('Adapt0:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('Adapt0:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('Adapt0:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('Adapt0:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('Adapt0:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('Adapt0:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('Adapt0:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('Adapt0:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('Adapt0:bruss', './ark_bruss.exe', p),
          ark.Test('Adapt0:heat1D', './ark_heat1D.exe', p),
          ark.Test('Adapt0:heat2D', './ark_heat2D.exe', p),
          ark.Test('Adapt0:hires', './ark_hires.exe', p),
          ark.Test('Adapt0:medakzo', './ark_medakzo.exe', p),
          ark.Test('Adapt0:orego', './ark_orego.exe', p),
          ark.Test('Adapt0:pollu', './ark_pollu.exe', p),
          ark.Test('Adapt0:ringmod', './ark_ringmod.exe', p),
          ark.Test('Adapt0:rober', './ark_rober.exe', p),
          ark.Test('Adapt0:robertson', './ark_robertson.exe', p),
          ark.Test('Adapt0:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('Adapt0:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('Adapt0:vdpol', './ark_vdpol.exe', p),
          ark.Test('Adapt0:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('Adapt0:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('Adapt0:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('Adapt0:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('Adapt0:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('Adapt0:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('Adapt0:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('Adapt0:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('Adapt0:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('Adapt0:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('Adapt0:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('Adapt0:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('Adapt0:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('Adapt0:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('Adapt0:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('Adapt0:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('Adapt0:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('Adapt0:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('Adapt0:tight-hires', './ark_hires.exe', p2),
          ark.Test('Adapt0:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('Adapt0:tight-orego', './ark_orego.exe', p2),
          ark.Test('Adapt0:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('Adapt0:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('Adapt0:tight-rober', './ark_rober.exe', p2),
          ark.Test('Adapt0:tight-robertson', './ark_robertson.exe', p2),
          ark.Test('Adapt0:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('Adapt0:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('Adapt0:tight-vdpol', './ark_vdpol.exe', p2),
          ark.Test('Adapt0:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### adaptivity method 1 tests ####
p = ark.SolParams(adapt_method=1)
tests = ( ark.Test('Adapt1:analytic', './ark_analytic.exe', p),
          ark.Test('Adapt1:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('Adapt1:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('Adapt1:system', './ark_analytic_sys.exe', p),
          ark.Test('Adapt1:brusselator', './ark_brusselator.exe', p),
          ark.Test('Adapt1:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('Adapt1:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('Adapt1:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('Adapt1:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('Adapt1:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('Adapt1:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('Adapt1:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('Adapt1:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('Adapt1:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('Adapt1:bruss', './ark_bruss.exe', p),
          ark.Test('Adapt1:heat1D', './ark_heat1D.exe', p),
          ark.Test('Adapt1:heat2D', './ark_heat2D.exe', p),
          ark.Test('Adapt1:hires', './ark_hires.exe', p),
          ark.Test('Adapt1:medakzo', './ark_medakzo.exe', p),
          ark.Test('Adapt1:pollu', './ark_pollu.exe', p),
          ark.Test('Adapt1:ringmod', './ark_ringmod.exe', p),
          ark.Test('Adapt1:rober', './ark_rober.exe', p),
          ark.Test('Adapt1:robertson', './ark_robertson.exe', p),
          ark.Test('Adapt1:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('Adapt1:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('Adapt1:vdpol', './ark_vdpol.exe', p),
          ark.Test('Adapt1:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('Adapt1:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('Adapt1:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('Adapt1:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('Adapt1:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('Adapt1:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('Adapt1:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('Adapt1:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('Adapt1:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('Adapt1:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('Adapt1:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('Adapt1:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('Adapt1:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('Adapt1:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('Adapt1:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('Adapt1:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('Adapt1:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('Adapt1:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('Adapt1:tight-hires', './ark_hires.exe', p2),
          ark.Test('Adapt1:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('Adapt1:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('Adapt1:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('Adapt1:tight-rober', './ark_rober.exe', p2),
          ark.Test('Adapt1:tight-robertson', './ark_robertson.exe', p2),
          ark.Test('Adapt1:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('Adapt1:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('Adapt1:tight-vdpol', './ark_vdpol.exe', p2),
          ark.Test('Adapt1:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### adaptivity method 2 tests ####
p = ark.SolParams(adapt_method=2)
tests = ( ark.Test('Adapt2:analytic', './ark_analytic.exe', p),
          ark.Test('Adapt2:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('Adapt2:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('Adapt2:system', './ark_analytic_sys.exe', p),
          ark.Test('Adapt2:brusselator', './ark_brusselator.exe', p),
          ark.Test('Adapt2:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('Adapt2:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('Adapt2:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('Adapt2:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('Adapt2:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('Adapt2:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('Adapt2:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('Adapt2:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('Adapt2:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('Adapt2:bruss', './ark_bruss.exe', p),
          ark.Test('Adapt2:heat1D', './ark_heat1D.exe', p),
          ark.Test('Adapt2:heat2D', './ark_heat2D.exe', p),
          ark.Test('Adapt2:hires', './ark_hires.exe', p),
          ark.Test('Adapt2:medakzo', './ark_medakzo.exe', p),
          ark.Test('Adapt2:orego', './ark_orego.exe', p),
          ark.Test('Adapt2:pollu', './ark_pollu.exe', p),
          ark.Test('Adapt2:ringmod', './ark_ringmod.exe', p),
          ark.Test('Adapt2:rober', './ark_rober.exe', p),
          ark.Test('Adapt2:robertson', './ark_robertson.exe', p),
          ark.Test('Adapt2:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('Adapt2:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('Adapt2:vdpol', './ark_vdpol.exe', p),
          ark.Test('Adapt2:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('Adapt2:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('Adapt2:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('Adapt2:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('Adapt2:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('Adapt2:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('Adapt2:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('Adapt2:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('Adapt2:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('Adapt2:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('Adapt2:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('Adapt2:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('Adapt2:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('Adapt2:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('Adapt2:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('Adapt2:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('Adapt2:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('Adapt2:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('Adapt2:tight-hires', './ark_hires.exe', p2),
          ark.Test('Adapt2:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('Adapt2:tight-orego', './ark_orego.exe', p2),
          ark.Test('Adapt2:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('Adapt2:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('Adapt2:tight-rober', './ark_rober.exe', p2),
          ark.Test('Adapt2:tight-robertson', './ark_robertson.exe', p2),
          ark.Test('Adapt2:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('Adapt2:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('Adapt2:tight-vdpol', './ark_vdpol.exe', p2),
          ark.Test('Adapt2:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### adaptivity method 3 tests ####
p = ark.SolParams(adapt_method=3)
tests = ( ark.Test('Adapt3:analytic', './ark_analytic.exe', p),
          ark.Test('Adapt3:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('Adapt3:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('Adapt3:system', './ark_analytic_sys.exe', p),
          ark.Test('Adapt3:brusselator', './ark_brusselator.exe', p),
          ark.Test('Adapt3:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('Adapt3:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('Adapt3:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('Adapt3:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('Adapt3:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('Adapt3:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('Adapt3:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('Adapt3:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('Adapt3:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('Adapt3:bruss', './ark_bruss.exe', p),
          ark.Test('Adapt3:heat1D', './ark_heat1D.exe', p),
          ark.Test('Adapt3:heat2D', './ark_heat2D.exe', p),
          ark.Test('Adapt3:hires', './ark_hires.exe', p),
          ark.Test('Adapt3:medakzo', './ark_medakzo.exe', p),
          ark.Test('Adapt3:orego', './ark_orego.exe', p),
          ark.Test('Adapt3:pollu', './ark_pollu.exe', p),
          ark.Test('Adapt3:ringmod', './ark_ringmod.exe', p),
          ark.Test('Adapt3:rober', './ark_rober.exe', p),
          ark.Test('Adapt3:robertson', './ark_robertson.exe', p),
          ark.Test('Adapt3:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('Adapt3:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('Adapt3:vdpol', './ark_vdpol.exe', p),
          ark.Test('Adapt3:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('Adapt3:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('Adapt3:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('Adapt3:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('Adapt3:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('Adapt3:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('Adapt3:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('Adapt3:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('Adapt3:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('Adapt3:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('Adapt3:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('Adapt3:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('Adapt3:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('Adapt3:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('Adapt3:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('Adapt3:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('Adapt3:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('Adapt3:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('Adapt3:tight-hires', './ark_hires.exe', p2),
          ark.Test('Adapt3:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('Adapt3:tight-orego', './ark_orego.exe', p2),
          ark.Test('Adapt3:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('Adapt3:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('Adapt3:tight-rober', './ark_rober.exe', p2),
          ark.Test('Adapt3:tight-robertson', './ark_robertson.exe', p2),
          ark.Test('Adapt3:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('Adapt3:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('Adapt3:tight-vdpol', './ark_vdpol.exe', p2),
          ark.Test('Adapt3:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### adaptivity method 4 tests ####
p = ark.SolParams(adapt_method=4)
tests = ( ark.Test('Adapt4:analytic', './ark_analytic.exe', p),
          ark.Test('Adapt4:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('Adapt4:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('Adapt4:system', './ark_analytic_sys.exe', p),
          ark.Test('Adapt4:brusselator', './ark_brusselator.exe', p),
          ark.Test('Adapt4:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('Adapt4:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('Adapt4:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('Adapt4:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('Adapt4:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('Adapt4:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('Adapt4:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('Adapt4:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('Adapt4:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('Adapt4:bruss', './ark_bruss.exe', p),
          ark.Test('Adapt4:heat1D', './ark_heat1D.exe', p),
          ark.Test('Adapt4:heat2D', './ark_heat2D.exe', p),
          ark.Test('Adapt4:hires', './ark_hires.exe', p),
          ark.Test('Adapt4:medakzo', './ark_medakzo.exe', p),
          ark.Test('Adapt4:orego', './ark_orego.exe', p),
          ark.Test('Adapt4:pollu', './ark_pollu.exe', p),
          ark.Test('Adapt4:ringmod', './ark_ringmod.exe', p),
          ark.Test('Adapt4:rober', './ark_rober.exe', p),
          ark.Test('Adapt4:robertson', './ark_robertson.exe', p),
          ark.Test('Adapt4:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('Adapt4:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('Adapt4:vdpol', './ark_vdpol.exe', p),
          ark.Test('Adapt4:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('Adapt4:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('Adapt4:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('Adapt4:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('Adapt4:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('Adapt4:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('Adapt4:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('Adapt4:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('Adapt4:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('Adapt4:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('Adapt4:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('Adapt4:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('Adapt4:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('Adapt4:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('Adapt4:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('Adapt4:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('Adapt4:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('Adapt4:tight-hires', './ark_hires.exe', p2),
          ark.Test('Adapt4:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('Adapt4:tight-orego', './ark_orego.exe', p2),
          ark.Test('Adapt4:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('Adapt4:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('Adapt4:tight-rober', './ark_rober.exe', p2),
          ark.Test('Adapt4:tight-robertson', './ark_robertson.exe', p2),
          ark.Test('Adapt4:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('Adapt4:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('Adapt4:tight-vdpol', './ark_vdpol.exe', p2) )
AllTests.append(tests)


#### adaptivity method 5 tests ####
p = ark.SolParams(adapt_method=5)
tests = ( ark.Test('Adapt5:analytic', './ark_analytic.exe', p),
          ark.Test('Adapt5:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('Adapt5:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('Adapt5:system', './ark_analytic_sys.exe', p),
          ark.Test('Adapt5:brusselator', './ark_brusselator.exe', p),
          ark.Test('Adapt5:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('Adapt5:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('Adapt5:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('Adapt5:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('Adapt5:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('Adapt5:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('Adapt5:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('Adapt5:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('Adapt5:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('Adapt5:bruss', './ark_bruss.exe', p),
          ark.Test('Adapt5:heat1D', './ark_heat1D.exe', p),
          ark.Test('Adapt5:heat2D', './ark_heat2D.exe', p),
          ark.Test('Adapt5:hires', './ark_hires.exe', p),
          ark.Test('Adapt5:medakzo', './ark_medakzo.exe', p),
          ark.Test('Adapt5:orego', './ark_orego.exe', p),
          ark.Test('Adapt5:pollu', './ark_pollu.exe', p),
          ark.Test('Adapt5:ringmod', './ark_ringmod.exe', p),
          ark.Test('Adapt5:rober', './ark_rober.exe', p),
          ark.Test('Adapt5:robertson', './ark_robertson.exe', p),
          ark.Test('Adapt5:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('Adapt5:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('Adapt5:vdpol', './ark_vdpol.exe', p),
          ark.Test('Adapt5:vdpolm', './ark_vdpolm.exe', p),
          ark.Test('Adapt5:tight-analytic', './ark_analytic.exe', p2),
          ark.Test('Adapt5:tight-nonlin', './ark_analytic_nonlin.exe', p2),
          ark.Test('Adapt5:tight-nonlin_back', './ark_analytic_nonlin_back.exe', p2),
          ark.Test('Adapt5:tight-system', './ark_analytic_sys.exe', p2),
          ark.Test('Adapt5:tight-brusselator', './ark_brusselator.exe', p2),
          ark.Test('Adapt5:tight-bruss_Ma', './ark_brusselator_Ma.exe', p2),
          ark.Test('Adapt5:tight-bruss_Mb', './ark_brusselator_Mb.exe', p2),
          ark.Test('Adapt5:tight-bruss_Mc', './ark_brusselator_Mc.exe', p2),
          ark.Test('Adapt5:tight-bruss_Md', './ark_brusselator_Md.exe', p2),
          ark.Test('Adapt5:tight-bruss_Me', './ark_brusselator_Me.exe', p2),
          ark.Test('Adapt5:tight-bruss1D', './ark_brusselator1D.exe', p2),
          ark.Test('Adapt5:tight-bruss1D_klu', './ark_brusselator1D_klu.exe', p2),
          ark.Test('Adapt5:tight-bruss1D_FEM', './ark_brusselator1D_FEM.exe', p2),
          ark.Test('Adapt5:tight-bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p2),
          ark.Test('Adapt5:tight-bruss', './ark_bruss.exe', p2),
          ark.Test('Adapt5:tight-heat1D', './ark_heat1D.exe', p2),
          ark.Test('Adapt5:tight-heat2D', './ark_heat2D.exe', p2),
          ark.Test('Adapt5:tight-hires', './ark_hires.exe', p2),
          ark.Test('Adapt5:tight-medakzo', './ark_medakzo.exe', p2),
          ark.Test('Adapt5:tight-orego', './ark_orego.exe', p2),
          ark.Test('Adapt5:tight-pollu', './ark_pollu.exe', p2),
          ark.Test('Adapt5:tight-ringmod', './ark_ringmod.exe', p2),
          ark.Test('Adapt5:tight-rober', './ark_rober.exe', p2),
          ark.Test('Adapt5:tight-robertson', './ark_robertson.exe', p2),
          ark.Test('Adapt5:tight-robertson_klu', './ark_robertson_klu.exe', p2),
          #ark.Test('Adapt5:tight-robertson_root', './ark_robertson_root.exe', p2),
          ark.Test('Adapt5:tight-vdpol', './ark_vdpol.exe', p2),
          ark.Test('Adapt5:tight-vdpolm', './ark_vdpolm.exe', p2) )
AllTests.append(tests)


#### predictor method 0 tests ####
p = ark.SolParams(predictor=0)
tests = ( ark.Test('Predict0:analytic', './ark_analytic.exe', p),
          ark.Test('Predict0:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('Predict0:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('Predict0:system', './ark_analytic_sys.exe', p),
          ark.Test('Predict0:brusselator', './ark_brusselator.exe', p),
          ark.Test('Predict0:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('Predict0:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('Predict0:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('Predict0:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('Predict0:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('Predict0:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('Predict0:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('Predict0:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('Predict0:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('Predict0:bruss', './ark_bruss.exe', p),
          ark.Test('Predict0:heat1D', './ark_heat1D.exe', p),
          ark.Test('Predict0:heat2D', './ark_heat2D.exe', p),
          ark.Test('Predict0:hires', './ark_hires.exe', p),
          ark.Test('Predict0:medakzo', './ark_medakzo.exe', p),
          ark.Test('Predict0:orego', './ark_orego.exe', p),
          ark.Test('Predict0:pollu', './ark_pollu.exe', p),
          ark.Test('Predict0:ringmod', './ark_ringmod.exe', p),
          ark.Test('Predict0:rober', './ark_rober.exe', p),
          ark.Test('Predict0:robertson', './ark_robertson.exe', p),
          ark.Test('Predict0:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('Predict0:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('Predict0:vdpol', './ark_vdpol.exe', p),
          ark.Test('Predict0:vdpolm', './ark_vdpolm.exe', p) )
AllTests.append(tests)


#### predictor method 1 tests ####
p = ark.SolParams(predictor=1)
tests = ( ark.Test('Predict1:analytic', './ark_analytic.exe', p),
          ark.Test('Predict1:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('Predict1:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('Predict1:system', './ark_analytic_sys.exe', p),
          ark.Test('Predict1:brusselator', './ark_brusselator.exe', p),
          ark.Test('Predict1:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('Predict1:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('Predict1:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('Predict1:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('Predict1:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('Predict1:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('Predict1:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('Predict1:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('Predict1:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('Predict1:heat1D', './ark_heat1D.exe', p),
          ark.Test('Predict1:heat2D', './ark_heat2D.exe', p),
          ark.Test('Predict1:hires', './ark_hires.exe', p),
          ark.Test('Predict1:medakzo', './ark_medakzo.exe', p),
          ark.Test('Predict1:rober', './ark_rober.exe', p),
          ark.Test('Predict1:robertson', './ark_robertson.exe', p),
          ark.Test('Predict1:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('Predict1:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('Predict1:vdpol', './ark_vdpol.exe', p),
          ark.Test('Predict1:vdpolm', './ark_vdpolm.exe', p) )
AllTests.append(tests)


#### predictor method 2 tests ####
p = ark.SolParams(predictor=2)
tests = ( ark.Test('Predict2:analytic', './ark_analytic.exe', p),
          ark.Test('Predict2:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('Predict2:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('Predict2:system', './ark_analytic_sys.exe', p),
          ark.Test('Predict2:brusselator', './ark_brusselator.exe', p),
          ark.Test('Predict2:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('Predict2:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('Predict2:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('Predict2:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('Predict2:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('Predict2:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('Predict2:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('Predict2:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('Predict2:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('Predict2:bruss', './ark_bruss.exe', p),
          ark.Test('Predict2:heat1D', './ark_heat1D.exe', p),
          ark.Test('Predict2:heat2D', './ark_heat2D.exe', p),
          ark.Test('Predict2:hires', './ark_hires.exe', p),
          ark.Test('Predict2:medakzo', './ark_medakzo.exe', p),
          ark.Test('Predict2:orego', './ark_orego.exe', p),
          ark.Test('Predict2:pollu', './ark_pollu.exe', p),
          ark.Test('Predict2:ringmod', './ark_ringmod.exe', p),
          ark.Test('Predict2:rober', './ark_rober.exe', p),
          ark.Test('Predict2:robertson', './ark_robertson.exe', p),
          ark.Test('Predict2:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('Predict2:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('Predict2:vdpol', './ark_vdpol.exe', p),
          ark.Test('Predict2:vdpolm', './ark_vdpolm.exe', p) )
AllTests.append(tests)


#### predictor method 3 tests ####
p = ark.SolParams(predictor=3)
tests = ( ark.Test('Predict3:analytic', './ark_analytic.exe', p),
          ark.Test('Predict3:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('Predict3:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('Predict3:system', './ark_analytic_sys.exe', p),
          ark.Test('Predict3:brusselator', './ark_brusselator.exe', p),
          ark.Test('Predict3:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('Predict3:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('Predict3:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('Predict3:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('Predict3:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('Predict3:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('Predict3:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('Predict3:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('Predict3:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('Predict3:bruss', './ark_bruss.exe', p),
          ark.Test('Predict3:heat1D', './ark_heat1D.exe', p),
          ark.Test('Predict3:heat2D', './ark_heat2D.exe', p),
          ark.Test('Predict3:hires', './ark_hires.exe', p),
          ark.Test('Predict3:medakzo', './ark_medakzo.exe', p),
          ark.Test('Predict3:pollu', './ark_pollu.exe', p),
          ark.Test('Predict3:ringmod', './ark_ringmod.exe', p),
          ark.Test('Predict3:rober', './ark_rober.exe', p),
          ark.Test('Predict3:robertson', './ark_robertson.exe', p),
          ark.Test('Predict3:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('Predict3:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('Predict3:vdpol', './ark_vdpol.exe', p),
          ark.Test('Predict3:vdpolm', './ark_vdpolm.exe', p) )
AllTests.append(tests)


#### predictor method 4 tests ####
p = ark.SolParams(predictor=4)
tests = ( ark.Test('Predict4:analytic', './ark_analytic.exe', p),
          ark.Test('Predict4:nonlin', './ark_analytic_nonlin.exe', p),
          ark.Test('Predict4:nonlin_back', './ark_analytic_nonlin_back.exe', p),
          ark.Test('Predict4:system', './ark_analytic_sys.exe', p),
          ark.Test('Predict4:brusselator', './ark_brusselator.exe', p),
          ark.Test('Predict4:bruss_Ma', './ark_brusselator_Ma.exe', p),
          ark.Test('Predict4:bruss_Mb', './ark_brusselator_Mb.exe', p),
          ark.Test('Predict4:bruss_Mc', './ark_brusselator_Mc.exe', p),
          ark.Test('Predict4:bruss_Md', './ark_brusselator_Md.exe', p),
          ark.Test('Predict4:bruss_Me', './ark_brusselator_Me.exe', p),
          ark.Test('Predict4:bruss1D', './ark_brusselator1D.exe', p),
          ark.Test('Predict4:bruss1D_klu', './ark_brusselator1D_klu.exe', p),
          ark.Test('Predict4:bruss1D_FEM', './ark_brusselator1D_FEM.exe', p),
          ark.Test('Predict4:bruss1D_FEM_klu', './ark_brusselator1D_FEM_klu.exe', p),
          ark.Test('Predict4:bruss', './ark_bruss.exe', p),
          ark.Test('Predict4:heat1D', './ark_heat1D.exe', p),
          ark.Test('Predict4:heat2D', './ark_heat2D.exe', p),
          ark.Test('Predict4:hires', './ark_hires.exe', p),
          ark.Test('Predict4:medakzo', './ark_medakzo.exe', p),
          ark.Test('Predict4:orego', './ark_orego.exe', p),
          ark.Test('Predict4:pollu', './ark_pollu.exe', p),
          ark.Test('Predict4:ringmod', './ark_ringmod.exe', p),
          ark.Test('Predict4:rober', './ark_rober.exe', p),
          ark.Test('Predict4:robertson', './ark_robertson.exe', p),
          ark.Test('Predict4:robertson_klu', './ark_robertson_klu.exe', p),
          #ark.Test('Predict4:robertson_root', './ark_robertson_root.exe', p),
          ark.Test('Predict4:vdpol', './ark_vdpol.exe', p),
          ark.Test('Predict4:vdpolm', './ark_vdpolm.exe', p) )
AllTests.append(tests)



#### flatten test list ####
all_tests = list(y for x in AllTests for y in x)



#### create test description file ####
AllTests_file = open('local_reg_tests.txt','w')
AllTests_file.write("# ARKode regression tests\n")
AllTests_file.write("#\n")
AllTests_file.write("# Each line indicates a separate test.  All lines have the format:\n")
AllTests_file.write("#   Test <name> <executable> <solver params>\n")
AllTests_file.write("# where \n")
AllTests_file.write("#   <name> is a string (no spaces allowed) describing the test problem\n")
AllTests_file.write("#   <executable> is a string containing the program to execute\n")
AllTests_file.write("#   <solver params> is a set of numbers representing a SolParams\n")
AllTests_file.write("#                   object (see arkode_tools.py)\n")
AllTests_file.write("#\n")
AllTests_file.write("# All non-test lines must begin with a \# sign.\n")
AllTests_file.write("#\n")
AllTests_file.write("# Daniel R. Reynolds\n")
AllTests_file.write("# SMU Mathematics\n")
AllTests_file.write("# June 2014\n")


#### run tests, storing results to dictionary ####
tstart = time()
numsuccess = 0
numfailed = 0
test_results = {};
for test in all_tests:

    # run test
    stats = test.Run(0); 

    # report on test success/failure
    if (stats.errfail == 1):
        numfailed += 1
        sys.stdout.write("  %30s \033[91m test failure \033[94m(runtime: %.2g s)\033[0m\n" 
                         % (test.name, stats.runtime))
    else:
        numsuccess += 1
        sys.stdout.write("  %30s \033[94m steps: %6i  runtime: %.2g s\033[0m\n" 
                         % (test.name, stats.nsteps, stats.runtime))

    # if test successful, store results to dictionary, and test descriptor to file
    if (stats.errfail == 0):
        test_results[test.name] = stats
        test.Write(AllTests_file)


# stop timer, and close test description file
tend = time()

#### close test description file, write standards dictionary ####
AllTests_file.close()
pickle.dump( test_results, open( "local_standards.dat", 'wb') )

#### output final statistics ####
sys.stdout.write("\nTest suite includes %i tests (%i failed/discarded); required %g seconds\n" 
                 % (numsuccess, numfailed, tend-tstart))

#### end of script ####
