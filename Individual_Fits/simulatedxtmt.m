function [xtf,mtf,xmtf] = simulatedxtmt(inj,xti,mti,V0,SyC)

mtf=mti;
xtf=xti;
Mcp=0;
for i=1:length(inj)
    mtfi=((V0-inj(i))*mti)/V0;
    if i==1
        injc=inj(i)*SyC/V0;
    else
        injc=((Mcp+inj(i)*SyC)-(inj(i)*(Mcp+inj(i)*SyC)/(V0+inj(i))))/V0;
    end
    
    mtf=vertcat(mtf,mtfi);
    xtf=vertcat(xtf,injc);
    
    mti=mtfi;
    Mcp=injc*V0;
end
xtf=xtf(2:length(xtf));
mtf=mtf(2:length(mtf));

xmtf=xtf./mtf;

end

