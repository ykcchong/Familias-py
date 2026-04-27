"""Internal port of the C++ Familias engine.

Module mapping (C++ -> Python):
    aldata.{h,cpp}, alsys.{h,cpp}     -> allele_system.py
    person.{h,cpp}, family.{h,cpp}    -> family.py
    odds.{h,cpp}, cutset.{h,cpp}      -> peeling.py
    Pedigree.{h,cpp}                  -> pedigree_core.py
    PedigreeList.{h,cpp}              -> pedigree_list.py
    pater.{h,cpp}, FamInterface.*     -> interface.py
    main.cpp                          -> (replaced by Python public API)
"""
from .interface import FamInterface

__all__ = ["FamInterface"]
