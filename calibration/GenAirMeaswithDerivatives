function [f,j] = GenAirMeaswithDerivatives(amp,a_PSG,b_PSG,a_PSA,b_PSA,PSG_theta,PSA_theta,PSA_LP,nSteps,LambdaList)
%a PSG is slope of PSG retardance 
% b PSG is slope intercept of PSG retardance
% same for PSA

for ii = 1:nLambda
for n = 1:64

 PSG_delta=a_PSG*LambdaList(n)+b_PSG;
 PSA_delta=a_PSA*LambdaList(n)+b_PSA;
    
a_12 = (cos(2*PSA_LP)*((cos(2*(n*5.625*4.91+PSA_theta)).^2) + (cos(PSA_delta)*sin(2*(n*5.625*4.91+PSA_theta)).^2))) + (sin(2*PSA_LP)*(((1-cos(PSA_delta))*cos(2*(n*5.625*4.91+PSA_theta))*sin(2*(n*5.625*4.91+PSA_theta)))));

a_13 = cos(2*PSA_LP)*((1-cos(PSA_delta))*cos(2*(n*5.625*4.91+PSA_theta))*sin(2*(n*5.625*4.91+PSA_theta)))+ ((sin(2*PSA_LP)*((cos(PSA_delta)*cos(2*(n*5.625*4.91+PSA_theta)).^2) + sin(2*(n*5.625*4.91+PSA_theta)).^2)));

a_14 = cos(2*PSA_LP)*(-sin(PSA_delta)*sin(2*(n*5.625*4.91+PSA_theta)))+(sin(2*PSA_LP)*(cos(2*(n*5.625*4.91+PSA_theta))*sin(PSA_delta)));

g_21 = cos(2*(n*5.625+PSG_theta)).^2 + (cos(PSG_delta)*(sin(2*(n*5.625+PSG_theta)).^2));

g_31 = (1-cos(PSG_delta))*cos(2*(n*5.625+PSG_theta))*sin(2*(n*5.625+PSG_theta));

g_41 = sin(PSG_delta)*sin(2*(n*5.625+PSG_theta));

%dpn_dPSG_delta = (a_12.*((cos(2*(n*5.625+PSG_theta)).^2)-(sin(PSG_delta)*sin(2*(n*5.625+PSG_theta)).^2)))+ (a_13.*(sin(PSG_delta)*cos(2*(n*5.625+PSG_theta))*sin(2*(n*5.625+PSG_theta))))+(a_14.*(cos(PSG_delta)*sin(2*(n*5.625+PSG_theta))));

ddelta_a_PSA = LambdaList(ii);
ddelta_a_PSG = LambdaList(ii);
ddelta_b_PSA = 1;
ddelta_b_PSG = 1;

dpn_dPSG_theta = (a_12.*((-2*sin(4*(n*5.625+PSG_theta)))+(cos(PSG_delta)*2*sin(4*(n*5.625+PSG_theta)))))+ (a_13.*((1-cos(PSG_delta))*(-2*sin(2*(n*5.625+PSG_theta)))*(2*cos(2*(n*5.625+PSG_theta)))))+(a_14.*(sin(PSG_delta)*2*cos(2*(n*5.625+PSG_theta)))); 

%dpn_dPSA_delta = ((((cos(2*PSA_LP))*((cos(2*(n*5.625*4.91+PSA_theta)).^2)+((-sin(PSA_delta))*(sin(2*(n*5.625*4.91+PSA_theta)).^2)))) + ((sin(2*PSA_LP)) * ((sin(PSA_delta))*(cos(2*(n*5.625*4.91+PSA_theta)))*(sin(2*(n*5.625*4.91+PSA_theta)))))).* g_21) + ((((cos(2*PSA_LP))*((1-sin(PSA_delta))*(cos(2*(n*5.625*4.91+PSA_theta)))*(sin(2*(n*5.625*4.91+PSA_theta))))) + ((sin(2*PSA_LP))*(((-sin(PSA_delta))*(cos(2*(n*5.625*4.91+PSA_theta)).^2)) + (sin(2*(n*5.625*4.91+PSA_theta).^2))))).*g_31) + ((((cos(2*(PSA_LP)))*((cos(PSA_delta))*(sin(2*(n*5.625*4.91+PSA_theta))))) + ((sin(2*PSA_LP))*((cos(2*(n*5.625*4.91+PSA_theta)))*(cos(PSA_delta))))).* g_41);

dpn_dPSA_theta = ((((cos(2*PSA_LP))*((-2*sin(4*(n*5.625*4.91+PSA_theta)))+((cos(PSA_delta))*(2*sin(4*(n*5.625*4.91+PSA_theta)))))) + ((sin(2*PSA_LP))*((1-cos(PSA_delta))*(-2*sin(2*(n*5.625*4.91+PSA_theta)))*(2*cos(2*(n*5.625*4.91+PSA_theta))))))* g_21) + ((((cos(2*PSA_LP))*((1-cos(PSA_delta))*(-2*sin(2*(n*5.625*4.91+PSA_theta)))*(2*cos(2*(n*5.625*4.91+PSA_theta))))) + ((sin(2*PSA_LP))*(((cos(PSA_delta))*(2*sin(4*(n*5.625*4.91+PSA_theta))))+(2*sin(4*(n*5.625*4.91+PSA_theta)))))).* g_31) + ((((cos(2*PSA_LP))*((-sin(PSA_delta))*(2*cos(2*(n*5.625*4.91+PSA_theta)))))+((sin(2*PSA_LP))*((-2*sin(2*(n*5.625*4.91+PSA_theta)))*(sin(PSA_delta))))).* g_41);

dpn_PSA_LP = ((-2*sin(2*PSA_LP)*((cos(2*(n*5.625*4.91+PSA_theta)).^2) +cos(PSA_delta)*sin(2*(n*5.625*4.91+PSA_theta)).^2)) + (2*cos(2*PSA_LP)*(((1-cos(PSA_delta))*cos(2*(n*5.625*4.91+PSA_theta))*sin(2*(n*5.625*4.91+PSA_theta)))) .* g_21)) + ((-2*sin(2*PSA_LP))*((1-cos(PSA_delta))*cos(2*(n*5.625*4.91+PSA_theta))*sin(2*(n*5.625*4.91+PSA_theta))) + (2*cos(2*PSA_LP)*((cos(PSA_delta)*cos(2*(n*5.625*4.91+PSA_theta)).^2) + sin(2*(n*5.625*4.91+PSA_theta)).^2)) .* g_31) + (((-2*sin(2*PSA_LP))*(-sin(PSA_delta)*sin(2*(n*5.625*4.91+PSA_theta)))) +((2*cos(2*PSA_LP))*(cos(2*(n*5.625*4.91+PSA_theta))*sin(PSA_delta))).* g_41);

% objective function values at x
f = GenAirMeasFromStepperParsV2(amp,a_PSG,b_PSG,a_PSA,b_PSA,PSG_theta,PSA_theta,PSA_LP,nSteps);
% Gradient of the objective function:
  if nargout  > 1
    j = [ddelta_a_PSG ddelta_b_PSG dpn_dPSG_theta ddelta_a_PSA ddelta_b_PSA dpn_dPSA_theta dpn_PSA_LP];
  end
end
end
end
