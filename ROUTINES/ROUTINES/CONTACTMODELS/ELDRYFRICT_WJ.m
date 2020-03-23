function [Fs, z, dFsdUs, dFsdUds, dFsdkxynmukth] = ELDRYFRICT_WJ(us, z, uds, kxynmukth, N0)
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
  dFsdkxynmukth = zeros(6, 9, Np);

				% 1. STICK (PREDICTION)
  Fs(1, :) = kxynmukth(1,:).*(us(1, :)-z(1, :));
  Fs(2, :) = kxynmukth(2,:).*(us(2, :)-z(2, :));
  Fs(3, :) = kxynmukth(3,:).*max(us(3, :)+N0./kxynmukth(3,:), 0);
  Fs(4, :) = kxynmukth(5, :).*us(4, :) + kxynmukth(6, :).*us(5, :);
  Fs(5, :) = kxynmukth(7, :).*us(4, :) + kxynmukth(8, :).*us(5, :);
  Fs(6, :) = kxynmukth(9, :).*us(6, :);
					    % Derivatives
  dFsdUs(1, 1, :) = kxynmukth(1,:); % fx,ux
  dFsdUs(2, 2, :) = kxynmukth(2,:); % fx,uy
  dFsdUs(3, 3, :) = kxynmukth(3,:); % fx,un
  dFsdUs(4, 4, :) = kxynmukth(5, :); % ftx,thx
  dFsdUs(4, 5, :) = kxynmukth(6, :); % ftx,thy
  dFsdUs(5, 4, :) = kxynmukth(7, :); % fty,thx
  dFsdUs(5, 5, :) = kxynmukth(8, :); % fty,thx
  dFsdUs(6, 6, :) = kxynmukth(9, :); % ftn,thn
					    % Parameter derivatives
  dFsdkxynmukth(1, 1, :) = us(1,:)-z(1,:); % fx,kx
  dFsdkxynmukth(2, 2, :) = us(2,:)-z(2,:); % fy,ky
  dFsdkxynmukth(3, 3, :) = us(3,:); % fn,kn
  dFsdkxynmukth(4, 5, :) = us(4, :); % ftx,ktx
  dFsdkxynmukth(4, 6, :) = us(5, :); % ftx,kty
  dFsdkxynmukth(5, 7, :) = us(4, :); % fty,ktx
  dFsdkxynmukth(5, 8, :) = us(5, :); % fty,kty
  dFsdkxynmukth(6, 9, :) = us(6, :); % ftn,ktn

				% 2. SEPARATION
  isep = find(Fs(3,:)==0);  % indices of separated points
  z(:, isep) = us(1:2, isep);
				% Everything is zero when separated
  Fs(:, isep) = 0;
  dFsdUs(:, :, isep) = 0;
  dFsdUds(:, :, isep) = 0;
  dFsdkxynmukth(:, :, isep) = 0;

				      % 3. SLIP
  fT = sqrt(sum(Fs(1:2, :).^2, 1));  % Tangential force magnitude
  fslip = kxynmukth(4,:).*Fs(3, :);        % Slip force magnitude

  islips = find(fT>fslip);       % indices of slipped points
  islips = setdiff(islips, isep);% Not interested in separated points - everything's zero there
  
  fT(fT<eps) = 1.0;  % Avoid dividing by zeros
  if ~isempty(islips)
				% Derivatives
    dFsdUs(1, 1, islips) = fslip(islips).*Fs(2,islips).^2.*kxynmukth(1,islips)./fT(islips).^3;% fx,ux
    dFsdUs(1, 2, islips) = -fslip(islips).*prod(Fs(1:2,islips),1).*kxynmukth(2,islips)./fT(islips).^3;% fx,uy
    dFsdUs(1, 3, islips) = kxynmukth(4,islips).*kxynmukth(3,islips).*Fs(1,islips)./fT(islips);% fx,un
    
    dFsdUs(2, 1, islips) = -fslip(islips).*prod(Fs(1:2,islips),1).*kxynmukth(1,islips)./fT(islips).^3;% fy,ux
    dFsdUs(2, 2, islips) = fslip(islips).*Fs(1,islips).^2.*kxynmukth(2,islips)./fT(islips).^3;% fy,uy
    dFsdUs(2, 3, islips) = kxynmukth(4,islips).*kxynmukth(3,islips).*Fs(2,islips)./fT(islips);% fy,un
    dFsdUs(4:end, :, islips) = 0;
				% Parameter Derivatives
    dFsdkxynmukth(1, 1, islips) = fslip(islips).*Fs(2,islips).^2.*(us(1,islips)-z(1,islips))./fT(islips).^3;% fx,kx
    dFsdkxynmukth(1, 2, islips) = -fslip(islips).*prod(Fs(1:2,islips),1).*(us(2,islips)-z(2,islips))./fT(islips).^3;% fx,ky
    dFsdkxynmukth(1, 3, islips) = kxynmukth(4,islips).*us(3,islips).*Fs(1,islips)./fT(islips);% fx,kn
    dFsdkxynmukth(1, 4, islips) = (fslip(islips)./kxynmukth(4,islips)).*Fs(1,islips)./fT(islips);% fx,mu

    dFsdkxynmukth(2, 1, islips) = -fslip(islips).*prod(Fs(1:2,islips),1).*(us(1,islips)-z(1,islips))./fT(islips).^3;% fy,kx  
    dFsdkxynmukth(2, 2, islips) = fslip(islips).*Fs(1,islips).^2.*(us(2,islips)-z(2,islips))./fT(islips).^3;% fy,ky
    dFsdkxynmukth(2, 3, islips) = kxynmukth(4,islips).*us(3,islips).*Fs(2,islips)./fT(islips);% fy,kn
    dFsdkxynmukth(2, 4, islips) = (fslip(islips)./kxynmukth(4,islips)).*Fs(2,islips)./fT(islips);% fy,mu
    dFsdkxynmukth(4:end, :, islips) = 0;

				% Update forces to lie on slip-cone
    Fs(1:2, islips) = Fs(1:2, islips).*repmat(fslip(islips)./fT(islips), 2, 1);
    Fs(4:end, islips) = 0;

    z(:, islips) = us(1:2, islips) - Fs(1:2, islips)./kxynmukth(1:2, islips);  % only x, y are hysteretic
  end
end
