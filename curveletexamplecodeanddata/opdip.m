classdef opdip < opSpot
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        m1,n1,dt,dx,p1,p2,p3,p4;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opdip(m1,n1,dt,dx,p1,p2,p3,p4)
            
            op = op@opSpot('opdip',m1*n1,m1*n1);
            op.cflag     = 0;
            op.linear    = 1;
            op.children  = [];
            op.sweepflag = true;
            op.dt         = dt;
            op.dx         = dx;
            op.p1         = p1;
            op.p2         = p2;
            op.p3         = p3;
            op.p4         = p4;
            op.m1         = m1;
            op.n1         = n1;
        end
    end
    
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = multiply(op,x,mode)
            if mode==1
                out = dipfilter(x,op.dt,op.dx,op.p1,op.p2,op.p3,op.p4,op.m1,op.n1);
            else
                out = dipfilter(x,op.dt,op.dx,op.p1,op.p2,op.p3,op.p4,op.m1,op.n1);
            end
        end %multiply
    end %protected methods
    
end %classdef

