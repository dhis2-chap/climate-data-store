from ecmwfapi import ECMWFDataServer

server = ECMWFDataServer()


server.retrieve({
    "class": "ti",
    "dataset": "tigge",
    "date": "2024-07-01/to/2024-07-29",
    "expver": "prod",
    "grid": "0.5/0.5",
    "levtype": "sfc",
    "origin": "ecmf",
    "param": "167",
    "step": "0/6/12/18/24/30/36/42/48/54/60/66/72/78/84/90/96/102/108/114/120/126/132/138/144/150/156/162/168/174/180/186/192/198/204/210/216/222/228/234/240/246/252/258/264/270/276/282/288/294/300/306/312/318/324/330/336/342/348/354/360",
    "time": "00:00:00",
    "type": "cf",
    "target": "output"
})
#server.retrieve({
#    "class": "rd",
#    "dataset": "research",
#    "date": "1989-01-01",
#    "expver": "gpfk",
#    "levelist": "10/50/100/200/300/400/500/700/850/925/1000",
#    "levtype": "pl",
#    "number": "0/1/2/3/4/5/6/7/8/9/10/11/12/13/14",
#    "param": "129.128/130.128/131.128/132.128/133.128/138.128/155.128",
#    "step": "96-264/264-432/432-600/600-768",
#    "stream": "enfo",
#    "target": "weekly_mean_plev_output.grib",
#    "time": "00:00:00",
#    "type": "fcmean"
#})