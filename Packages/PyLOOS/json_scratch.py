import json
test_string = """{
    "noe": {
        "atom": "one",
        "atom": "two",
        "volume": [0,1,2]
    },
    "noe" : {
        "atom" : "one",
        "atom" : "two",
        "volume" : [0,2,3]
        },
    "noe" : {
        "atom" : "two",
        "atom" : "three",
        "volume" : [2,4,6]
        },
    "noe" : {
        "atom" : "two",
        "atom" : "four",
        "volume" : [2,3,6]
        },
    "noe" : {
        "atom" : "one",
        "atom" : "four",
        "volume" : [3, 6, 9]
        }
}
"""

obj = json.loads(test_string)
print(obj)
