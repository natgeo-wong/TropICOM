function solveanalytic(
    model :: SimpleIslandModel
)

    atm  = model.atm
    sfc  = model.sfc
    
    return calculateanalytic(model,atm,sfc)

end