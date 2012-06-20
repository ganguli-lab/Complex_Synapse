numfail=0;
t=0:0.2:10;
x=0:1;
for i=1:100
    [s,env]=CheatCurves(20,t,x,0.3,'cheatq',false,'cheatW',true);
    if any(s(end,2:end)>env(2:end))
        plot(t,real([s;env]))
        legend([cellstr(num2str(x')); {'envelope'}]);
        pause;
        numfail=numfail+1;
    end%if
    if mod(i,10)==0
        disp(i);
    end%if
end%for i
disp(['number of failures = ' int2str(numfail)]);