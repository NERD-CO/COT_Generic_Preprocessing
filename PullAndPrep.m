% Run all

function PullAndPrep(sstrs)
for i = 1:length(sstrs)
    sstr = sstrs{i};

    PullFromSynology(sstr);
    AutoGenKey(sstr);
    disp('Converting .mat AO files to .plx for sorting');
    Batch_AO2plx(sstr);
    GenRunDLCPy(sstr);
    disp('Good job. Now go sort those units and DL that C!');
end