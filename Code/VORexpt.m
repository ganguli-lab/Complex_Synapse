function [ ph ] = VORexpt( Wp,Wm,w,t,df,tchange,varargin )
%ph=VOREXPT(Wp,Wm,w,t,df,tchange,...) Simulation of VOR expt
%   T = time values
%   WP = potentiation transition rates
%   WM = depression transition rates
%   DF = Change in fraction of potentiation transitions
%   w = Weights of states (+/-1)

error(CheckSize(Wp,@ismat));%matrix
error(CheckSize(Wp,@issquare));%square
error(CheckSize(Wm,@(x)samesize(Wp,x),'samesize(Wp)'));%also square matrix of same size
error(CheckSize(df,@isscalar));
error(CheckValue(df,@(x) inrange(x,-1,1),'inrange(0,1)'));%fp in [0,1]
error(CheckSize(w,@iscol));
error(CheckValue(w,@(x) all(x.^2==1),'all w = +/-1'));

fp=0.5;
LinSpec='b';
Parent=gca;

varargin=assignApplicable(varargin);

S=LearningCurve(fp*Wp+(1-fp)*Wm,w,t,(fp+df)*Wp+(1-fp-df)*Wm);
ph(1)=plot(t,S(1)-S,LinSpec,'Parent',Parent,varargin{:});
hold(Parent,'on');
S=LearningCurve(fp*Wp+(1-fp)*Wm,w,t,(fp-df)*Wp+(1-fp+df)*Wm,tchange,(fp+df)*Wp+(1-fp-df)*Wm);
ph(2)=plot(t,S(1)-S,LinSpec,'Parent',Parent,varargin{:});

end

