    network{zz}=-0.939998046792114*ones(scale{zz}); %#ok<AGROW> % depolarization array (starting state)
    fired{zz}=zeros(scale{zz});  %#ok<AGROW> % cells firing this iteration
    %refract{zz}=rand(scale{zz}); %#ok<AGROW> % refractory firing resistance
    potency{zz}=zeros(scale{zz}); %#ok<AGROW> % neuro transmiter suply
    input{zz}=zeros(scale{zz}); %#ok<AGROW> %
    umm{zz}=-2.821443297088526*ones(scale{zz}); %#ok<AGROW>
    u{zz}=-2.8214*ones(scale{zz}); %#ok<AGROW>
    cond1=zeros(scale{zz});
    cond2=zeros(scale{zz});%#ok<AGROW>

    cond3=zeros(scale{zz});
