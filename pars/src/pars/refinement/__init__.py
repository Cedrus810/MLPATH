"""Refinement of selected representatives (design sections 18 and 43 MVP-5).

The refinement interface is backend-independent: PARS defines the contract,
backends implement it. ASE's NEB is the first adapter; string method and
saddle refinement can be added without touching the selection pipeline.
"""
