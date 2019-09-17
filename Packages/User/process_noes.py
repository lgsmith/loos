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
    def __init__(self, noe_list: List[Noe], pairs: dict = {}, meta: List[str] = []):
        self.list = noe_list
        self.pairs = pairs
        self.meta = meta

    @classmethod
    def from_json_dict(cls, json_obj: list) -> "NoeExperiment":
        return cls(
            list=[Noe.from_json_dict(n) for n in json_obj["list"]],
            pairs=json_obj["pairs"],
            meta=json_obj["meta"],
        )

    @classmethod
    def from_json_file(cls, filename: str) -> "NoeExperiment":
        with open(filename, "r") as f:
            expt_json = json.load(f)
            return cls.from_json_dict(expt_json)

    @classmethod
    def from_loos(cls, filename: str) -> "NoeExperiment":
        noe_list = []
        meta = []
        with open(filename, "r") as f:
            for line in f:
                l = line.split()
                if l[0] == "#":
                    meta.append(l)
                else:
                    pair = (
                        dict(resname=l[2], resid=int(l[3]), name=l[4], index=int(l[5])),
                        dict(resname=l[6], resid=int(l[7]), name=l[8], index=int(l[9])),
                    )
                    noe_list.append(Noe(float(l[0]), float(l[1]), pair))
        return cls(noe_list=noe_list, meta=meta)

    def lookup_pairs(self, key_pattern=["resname", "resid", "name"]):
        for noe in self.list:
            key = tuple(at[k] for at in noe.pair for k in key_pattern)
            self.pairs[key] = noe

    # save this NoeExperiment to this filename
    def save(self, filename):
        # a deep copy, I believe?
        dump_dict = {}
        dump_dict.update(self.__dict__)
        dump_dict["list"] = [n.__dict__ for n in self.list]
        with open(filename, "w") as f:
            json.dump(dump_dict, f)


# provide a simple parser for loos files, partly as smoke test
if __name__ == "__main__":
    from sys import argv

    expt = NoeExperiment.from_loos(argv[1])
    expt.save(argv[2])
