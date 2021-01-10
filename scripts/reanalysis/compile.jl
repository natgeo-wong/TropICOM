using DrWatson
@quickactivate "TroPrecLS"

include(srcdir("reanalysis.jl"))
include(srcdir("diurnal.jl"))

compilesavesfc.([
    "ssr","str","sshf","slhf","tsr","ttr",
    "hcc","mcc","lcc","tcc","sst","skt","t2m",
    "tcw","tcwv","lsm","sp"
])

compilesavesfchour.([
    "ssr","str","sshf","slhf","tsr","ttr",
    "hcc","mcc","lcc","tcc","sst","skt","t2m",
    "tcw","tcwv","sp"
])

compileresavediurnalsfc.([
    "ssr","str","sshf","slhf","tsr","ttr",
    "hcc","mcc","lcc","tcc","sst","skt","t2m",
    "tcw","tcwv","sp"
])

compilesavepre.(["t_air","w_air","z_air","q_air","r_air","cc_air","clwc_air","ciwc_air"])

compilesaveprehour("cc_air",levels=[
    50,70,100,125,150,175,200,225,
    250,300,350,400,450,500,550,600,650,700,750,
    775,800,825,850,875,900,925,950,975,1000
]); compileresavediurnalpre("cc_air")

compilesaveprehour("clwc_air",levels=[
    50,70,100,125,150,175,200,225,
    250,300,350,400,450,500,550,600,650,700,750,
    775,800,825,850,875,900,925,950,975,1000
]); compileresavediurnalpre("cc_air")

compilesaveprehour("clwc_air",levels=[
    50,70,100,125,150,175,200,225,
    250,300,350,400,450,500,550,600,650,700,750,
    775,800,825,850,875,900,925,950,975,1000
]); compileresavediurnalpre("clwc_air")

compilesaveprehour("q_air",levels=[
    50,70,100,125,150,175,200,225,
    250,300,350,400,450,500,550,600,650,700,750,
    775,800,825,850,875,900,925,950,975,1000
]); compileresavediurnalpre("q_air")

compilesaveprehour("t_air",levels=[
    50,70,100,125,150,175,200,225,
    250,300,350,400,450,500,550,600,650,700,750,
    775,800,825,850,875,900,925,950,975,1000
]); compileresavediurnalpre("t_air")

compilesavesfcfeb.([
    "ssr","str","sshf","slhf"
])
