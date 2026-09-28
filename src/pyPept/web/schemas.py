"""Request limits and notation choices shared by the CABILN handlers."""

from __future__ import annotations

from typing import Annotated, Literal

from pydantic import BaseModel, Field, StringConstraints

Notation = Annotated[
    str, StringConstraints(strip_whitespace=True, min_length=1, max_length=20000)
]
Symbol = Annotated[
    str,
    StringConstraints(
        strip_whitespace=True, min_length=1, max_length=100, pattern=r"^[A-Za-z0-9_]+$"
    ),
]
Slot = Annotated[int, Field(ge=1, le=64)]
Label = Annotated[str, StringConstraints(min_length=1, max_length=100)]


class _ImageReq(BaseModel):
    width: int = Field(default=960, ge=64, le=4096)
    height: int = Field(default=680, ge=64, le=4096)


class _CabilnReq(_ImageReq):
    cabiln: Notation
    seed: int = Field(default=0, ge=0, le=2**31 - 1)


class _SmilesReq(_ImageReq):
    smiles: Notation


class _VerifyReq(BaseModel):
    smiles: Notation
    cabiln: Notation


class _PreviewReq(_ImageReq):
    smiles: Notation
    width: int = Field(default=640, ge=64, le=4096)
    height: int = Field(default=440, ge=64, le=4096)


class _RegisterReq(BaseModel):
    chuckles: Notation
    chem_types: dict[Slot, Label]
    leaving: dict[Slot, Label]
    abbr: Symbol
    name: Annotated[
        str, StringConstraints(strip_whitespace=True, min_length=1, max_length=200)
    ]
    type: Literal["aa", "cap", "linker"] = "aa"
    subtype: Label = "modified"


class _ConvertReq(BaseModel):
    cabiln: Notation
    target: Literal["bracket", "branch"]


class _InsertBondReq(BaseModel):
    cabiln: str = Field(max_length=20000)
    host_residue_idx: int = Field(ge=0)
    new_abbr: Symbol
    r_host: Slot
    r_new: Slot
    target_residue_idx: int = Field(default=-1, ge=-1)


class _InsertBackboneReq(BaseModel):
    cabiln: str = Field(max_length=20000)
    after_idx: int = Field(ge=0)
    new_abbr: Symbol


class _SmilesToCabilnReq(BaseModel):
    smiles: Notation
    notation: Literal["percent", "bracket"] = "percent"


class _ToCabilnReq(BaseModel):
    input: Notation
    notation: Literal["percent", "bracket"] = "percent"


class _ValidateBondReq(BaseModel):
    chem_type_a: Label
    chem_type_b: Label
    abbr_a: str = Field(default="", max_length=100)
    slot_a: int = Field(default=0, ge=0, le=64)
    abbr_b: str = Field(default="", max_length=100)
    slot_b: int = Field(default=0, ge=0, le=64)


class _ReferenceReq(_ImageReq):
    input: Notation


class _MolBlockReq(_ImageReq):
    mol_block: str = Field(min_length=1, max_length=1000000)
