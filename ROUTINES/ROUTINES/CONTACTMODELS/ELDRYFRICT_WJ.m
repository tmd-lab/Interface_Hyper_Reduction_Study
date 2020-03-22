function [Fs, z, dFsdUs, dFsdUds, dFsdkxynmu] = ELDRYFRICT_WJ(us, z, uds, kxynmu, I2, N0)
%ELDRYFRICT_WJ returns the forces and jacobians of elastic dry
%friction for a whole joint formulation
%
% USAGE:
% ------
%
% INPUTS:
% -------
%
% OUTPUTS:
% --------
  Np = size(us, 2);
  Fs = zeros(6, Np);
  dFsdUs = zeros(6, 6, Np);
  dFsdUds = zeros(6, 6, Np);
  dFsdkxynmu = zeros(6, 4, Np);

				% 1. STICK (PREDICTION)
  Fs(1, :) = kxynmu(1,:).*(us(1, :)-z(1, :));
  Fs(2, :) = kxynmu(2,:).*(us(2, :)-z(2, :));
  Fs(3, :) = kxynmu(3,:).*max(us(3, :)+N0./kxynmu(3,:), 0);
  Fs(4, :) = I2(1, :).*kxynmu(1, :).*us(4, :) + I2(2, :).*kxynmu(2, :).*us(5, :);
  Fs(5, :) = I2(2, :).*kxynmu(1, :).*us(4, :) + I2(3, :).*kxynmu(2, :).*us(5, :);
  Fs(6, :) = I2(4, :).*kxynmu(3, :).*us(6, :);
					    % Derivatives
  dFsdUs(1, 1, :) = kxynmu(1,:); % fx,ux
  dFsdUs(2, 2, :) = kxynmu(2,:); % fx,uy
  dFsdUs(3, 3, :) = kxynmu(3,:); % fx,un
  dFsdUs(4, 1, :) = I2(1, :).*kxynmu(1, :); % ftx,thx
  dFsdUs(4, 2, :) = I2(2, :).*kxynmu(2, :); % ftx,thy
  dFsdUs(5, 1, :) = I2(2, :).*kxynmu(1, :); % fty,thx
  dFsdUs(5, 2, :) = I2(3, :).*kxynmu(2, :); % fty,thx
  dFsdUs(6, 6, :) = I2(4, :).*kxynmu(3, :); % ftn,thn
					    % Parameter derivatives
  dFsdkxynmu(1, 1, :) = us(1,:)-z(1,:); % fx,kx
  dFsdkxynmu(2, 2, :) = us(2,:)-z(2,:); % fy,ky
  dFsdkxynmu(3, 3, :) = us(3,:); % fn,kn
  dFsdxynmu(4, 1, :) = I2(1, :).*us(4, :); % ftx,kx
  dFsdxynmu(4, 2, :) = I2(2, :).*us(5, :); % ftx,ky
  dFsdxynmu(5, 1, :) = I2(2, :).*us(4, :); % fty,kx
  dFsdxynmu(5, 2, :) = I2(3, :).*us(5, :); % fty,ky
  dFsdxynmu(6, 3, :) = I2(4, :).*us(6, :); % ftn,kn

				% 2. SEPARATION
  isep = find(Fs(3,:)==0);  % indices of separated points
  z(:, isep) = us(1:2, isep);
				% Everything is zero when separated
  Fs(:, isep) = 0;
  dFsdUs(:, :, isep) = 0;
  dFsdUds(:, :, isep) = 0;
  dFsdkxynmu(:, :, isep) = 0;

				      % 3. SLIP
  fT = sqrt(sum(Fs(1:2, :).^2, 1));  % Tangential force magnitude
  fslip = kxynmu(4,:).*Fs(3, :);        % Slip force magnitude

  islips = find(fT>fslip);       % indices of slipped points
  islips = setdiff(islips, isep);% Not interested in separated points - everything's zero there
  
  fT(fT<eps) = 1.0;  % Avoid dividing by zeros
  if length(islips) ~= 0
				% Derivatives
    dFsdUs(1, 1, islips) = fslip(islips).*Fs(2,islips).^2.*kxynmu(1,islips)./fT(islips).^3;% fx,ux
    dFsdUs(1, 2, islips) = -fslip(islips).*prod(Fs(1:2,islips),1).*kxynmu(2,islips)./fT(islips).^3;% fx,uy
    dFsdUs(1, 3, islips) = kxynmu(4,islips).*kxynmu(3,islips).*Fs(1,islips)./fT(islips);% fx,un
    
    dFsdUs(2, 1, islips) = -fslip(islips).*prod(Fs(1:2,islips),1).*kxynmu(1,islips)./fT(islips).^3;% fy,ux
    dFsdUs(2, 2, islips) = fslip(islips).*Fs(1,islips).^2.*kxynmu(2,islips)./fT(islips).^3;% fy,uy
    dFsdUs(2, 3, islips) = kxynmu(4,islips).*kxynmu(3,islips).*Fs(2,islips)./fT(islips);% fy,un
    dFsdUs(4:end, :, islips) = 0;
				% Parameter Derivatives
    dFsdkxynmu(1, 1, islips) = fslip(islips).*Fs(2,islips).^2.*(us(1,islips)-z(1,islips))./fT(islips).^3;% fx,kx
    dFsdkxynmu(1, 2, islips) = -fslip(islips).*prod(Fs(1:2,islips),1).*(us(2,islips)-z(2,islips))./fT(islips).^3;% fx,ky
    dFsdkxynmu(1, 3, islips) = kxynmu(4,islips).*us(3,islips).*Fs(1,islips)./fT(islips);% fx,kn
    dFsdkxynmu(1, 4, islips) = (fslip(islips)./kxynmu(4,islips)).*Fs(1,islips)./fT(islips);% fx,mu

    dFsdkxynmu(2, 1, islips) = -fslip(islips).*prod(Fs(1:2,islips),1).*(us(1,islips)-z(1,islips))./fT(islips).^3;% fy,kx  
    dFsdkxynmu(2, 2, islips) = fslip(islips).*Fs(1,islips).^2.*(us(2,islips)-z(2,islips))./fT(islips).^3;% fy,ky
    dFsdkxynmu(2, 3, islips) = kxynmu(4,islips).*us(3,islips).*Fs(2,islips)./fT(islips);% fy,kn
    dFsdkxynmu(2, 4, islips) = (fslip(islips)./kxynmu(4,islips)).*Fs(2,islips)./fT(islips);% fy,mu
    dFsdkxynmu(4:end, :, islips) = 0;

				% Update forces to lie on slip-cone
    Fs(1:2, islips) = Fs(1:2, islips).*repmat(fslip(islips)./fT(islips), 2, 1);
    Fs(4:end, islips) = 0;

    z(:, islips) = us(1:2, islips) - Fs(1:2, islips)./kxynmu(1:2, islips);  % only x, y are hysteretic
  end
end
