function Q = q_unpack(q,bndy,N_s,N_c)

    % Redfine boundary:
    Q = nan(N_s,N_c);
    Q(~bndy) = q;

    % Define Left BC (which is actually the top row of the matrix):
%     Q(1,:) = Q(2,:);
    Q(1,:) = 4/3*Q(2,:)-1/3*Q(3,:);

    % Define Right BC (which is the bottom row of the matrix):
    Q(end,:) = 1;
    
end
  
  