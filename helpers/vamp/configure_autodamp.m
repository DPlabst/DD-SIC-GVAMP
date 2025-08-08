function autoDampCFG = configure_autodamp(nvArgs, ind_sic)

    if nvArgs.aD == 1
        autoDampCFG.maxBadSteps = 3;
        autoDampCFG.maxstepDecr = 0.3;
        autoDampCFG.dampoffs = 0;
        autoDampCFG.stepWindow = max(nvArgs.dampcfgStepW(1), ...
            round(nvArgs.dampcfgStepW(2) * exp(- nvArgs.dampcfgStepW(3) * (ind_sic - 1))));
        autoDampCFG.stepMin = 0.01;
        autoDampCFG.stepMax = 1;
        autoDampCFG.stepTol = 0.01;
        autoDampCFG.stepInc = 1.1;
        autoDampCFG.stepDec = 0.5;
    else
        autoDampCFG.maxBadSteps = inf;
        autoDampCFG.maxstepDecr = 1;
        autoDampCFG.dampoffs = inf;
        autoDampCFG.stepWindow = inf;
        autoDampCFG.stepMin = 0.05;
        autoDampCFG.stepMax = 1;
        autoDampCFG.stepTol = 0.05;
        autoDampCFG.stepInc = 1;
        autoDampCFG.stepDec = 1;
    end

end
