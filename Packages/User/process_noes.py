import json
from typing import Tuple, List


class Noe:
    def __init__(self, volume: float, mean_dist: float, pair: Tuple[dict, dict]):
        self.volume = volume
        self.mean_dist = mean_dist
        self.pair = pair

    # taken from:
    # https://stavshamir.github.io/python/2018/05/26/overloading-constructors-in-python.html
    @classmethod
    def from_json_dict(cls, noe_dict: dict) -> "Noe":
        return cls(
            volume=noe_dict["volume"],
            mean_dist=noe_dict["mean_dist"],
            pair=tuple(noe_dict["pair"]),
        )


class NoeExperiment:
    def __init__(self, noe_list: List[Noe]):
        self.list = noe_list
        self.pairs = {}

    @classmethod
    def from_json_list(cls, json_obj: list) -> "NoeExperiment":
        return cls([Noe.from_json_dict(n) for n in json_obj])

    @classmethod
    def from_json_file(cls, filename: str) -> "NoeExperiment":
        with open(filename, "r") as f:
            expt_json = json.load(f)
            return cls.from_json_list(expt_json)

    def lookup_pairs(self, key_pattern=["resname", "resid", "name"]):
        for noe in self.list:
            key = tuple(at[k] for at in noe.pair for k in key_pattern)
            self.pairs[key] = noe

    # save this NoeExperiment to this filename
    def save(self, filename):
        with open(filename, 'w') as f:
            json.dump(f, [n.__dict__ for n in self.list])