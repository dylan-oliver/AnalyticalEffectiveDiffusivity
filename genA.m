function A = genA(grid,right_interface_A0,right_interface_An,right_interface_Bn,top_interface_A0,top_interface_An,top_interface_Bn)

L = grid.L;
W = grid.W;

Nx = grid.Nx;
Ny = grid.Ny;

index_right = grid.index_right;
index_top = grid.index_top;

total_entries = (2*(L-2)*W + 2*(L-1)*W + (2*(L-2) + 2*(L-1))*W*Nx + 4*(L-1)*(W-1)*Ny)*(Nx+1) ...
    + (4*(L-1)*(W-1) + 4*(L-1)*(W-1)*Nx + (2*(W-2) + 2*(W-1))*L*Ny)*Ny;

rows = zeros(total_entries,1);
cols = rows;
vals = rows;

entry = 0;

% Right-side boundary continuity conditions
row_index = @(i,j) (j - 1 + (i-1)*(L-1))*(Nx+1);
inc = 1:(Nx+1);

% Bottom-left corner block
i = 1;
j = 1;

row_entries = row_index(i,j) + inc;

% B0(i,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j);
vals(entries) = right_interface_A0(i,j,2);
entry = entry + Nx + 1;

% A0(i,j+1) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j);
vals(entries) = right_interface_A0(i,j,3);
entry = entry + Nx + 1;

% B0(i,j+1) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j+1);
vals(entries) = right_interface_A0(i,j,4);
entry = entry + Nx + 1;

for k = 1:Nx

    % Bn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j) + k;
    vals(entries) = right_interface_An(k,:,i,j,2);
    entry = entry + Nx + 1;

    % An(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j) + k;
    vals(entries) = right_interface_An(k,:,i,j,3);
    entry = entry + Nx + 1;

    % Bn(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j+1) + k;
    vals(entries) = right_interface_An(k,:,i,j,4);
    entry = entry + Nx + 1;

end

for k = 1:Ny

    % Dn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j) + k;
    vals(entries) = right_interface_Bn(k,:,i,j,2);
    entry = entry + Nx + 1;

    % Dn(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j+1) + k;
    vals(entries) = right_interface_Bn(k,:,i,j,4);
    entry = entry + Nx + 1;

end


% Bottom-right corner block
i = 1;
j = L-1;

row_entries = row_index(i,j) + inc;

% A0(i,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j-1);
vals(entries) = right_interface_A0(i,j,1);
entry = entry + Nx + 1;

% B0(i,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j);
vals(entries) = right_interface_A0(i,j,2);
entry = entry + Nx + 1;

% A0(i,j+1) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j);
vals(entries) = right_interface_A0(i,j,3);
entry = entry + Nx + 1;

for k = 1:Nx

    % An(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j-1) + k;
    vals(entries) = right_interface_An(k,:,i,j,1);
    entry = entry + Nx + 1;

    % Bn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j) + k;
    vals(entries) = right_interface_An(k,:,i,j,2);
    entry = entry + Nx + 1;

    % An(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j) + k;
    vals(entries) = right_interface_An(k,:,i,j,3);
    entry = entry + Nx + 1;

end

for k = 1:Ny

    % Dn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j) + k;
    vals(entries) = right_interface_Bn(k,:,i,j,2);
    entry = entry + Nx + 1;

    % Dn(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j+1) + k;
    vals(entries) = right_interface_Bn(k,:,i,j,4);
    entry = entry + Nx + 1;

end


% Top-left corner block
i = W;
j = 1;

row_entries = row_index(i,j) + inc;

% B0(i,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j);
vals(entries) = right_interface_A0(i,j,2);
entry = entry + Nx + 1;

% A0(i,j+1) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j);
vals(entries) = right_interface_A0(i,j,3);
entry = entry + Nx + 1;

% B0(i,j+1) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j+1);
vals(entries) = right_interface_A0(i,j,4);
entry = entry + Nx + 1;

for k = 1:Nx

    % Bn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j) + k;
    vals(entries) = right_interface_An(k,:,i,j,2);
    entry = entry + Nx + 1;

    % An(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j) + k;
    vals(entries) = right_interface_An(k,:,i,j,3);
    entry = entry + Nx + 1;

    % Bn(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j+1) + k;
    vals(entries) = right_interface_An(k,:,i,j,4);
    entry = entry + Nx + 1;

end

for k = 1:Ny

    % Cn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i-1,j) + k;
    vals(entries) = right_interface_Bn(k,:,i,j,1);
    entry = entry + Nx + 1;

    % Cn(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i-1,j+1) + k;
    vals(entries) = right_interface_Bn(k,:,i,j,3);
    entry = entry + Nx + 1;

end


% Top-right corner block
i = W;
j = L-1;

row_entries = row_index(i,j) + inc;

% A0(i,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j-1);
vals(entries) = right_interface_A0(i,j,1);
entry = entry + Nx + 1;

% B0(i,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j);
vals(entries) = right_interface_A0(i,j,2);
entry = entry + Nx + 1;

% A0(i,j+1) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j);
vals(entries) = right_interface_A0(i,j,3);
entry = entry + Nx + 1;

for k = 1:Nx

    % An(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j-1) + k;
    vals(entries) = right_interface_An(k,:,i,j,1);
    entry = entry + Nx + 1;

    % Bn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j) + k;
    vals(entries) = right_interface_An(k,:,i,j,2);
    entry = entry + Nx + 1;

    % An(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j) + k;
    vals(entries) = right_interface_An(k,:,i,j,3);
    entry = entry + Nx + 1;

end

for k = 1:Ny

    % Cn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i-1,j) + k;
    vals(entries) = right_interface_Bn(k,:,i,j,1);
    entry = entry + Nx + 1;

    % Cn(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i-1,j+1) + k;
    vals(entries) = right_interface_Bn(k,:,i,j,3);
    entry = entry + Nx + 1;

end


for j = 2:L-2
    % Bottom side blocks (not including corner blocks)
    i = 1;

    row_entries = row_index(i,j) + inc;

    % A0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j-1);
    vals(entries) = right_interface_A0(i,j,1);
    entry = entry + Nx + 1;

    % B0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j);
    vals(entries) = right_interface_A0(i,j,2);
    entry = entry + Nx + 1;

    % A0(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j);
    vals(entries) = right_interface_A0(i,j,3);
    entry = entry + Nx + 1;

    % B0(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j+1);
    vals(entries) = right_interface_A0(i,j,4);
    entry = entry + Nx + 1;

    for k = 1:Nx
        % An(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j-1) + k;
        vals(entries) = right_interface_An(k,:,i,j,1);
        entry = entry + Nx + 1;

        % Bn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j) + k;
        vals(entries) = right_interface_An(k,:,i,j,2);
        entry = entry + Nx + 1;

        % An(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j) + k;
        vals(entries) = right_interface_An(k,:,i,j,3);
        entry = entry + Nx + 1;

        % Bn(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j+1) + k;
        vals(entries) = right_interface_An(k,:,i,j,4);
        entry = entry + Nx + 1;
    end

    for k = 1:Ny
        % Dn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,2);
        entry = entry + Nx + 1;

        % Dn(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j+1) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,4);
        entry = entry + Nx + 1;
    end


    % Top side blocks (not including corner blocks)
    i = W;

    row_entries = row_index(i,j) + inc;

    % A0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j-1);
    vals(entries) = right_interface_A0(i,j,1);
    entry = entry + Nx + 1;

    % B0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j);
    vals(entries) = right_interface_A0(i,j,2);
    entry = entry + Nx + 1;

    % A0(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j);
    vals(entries) = right_interface_A0(i,j,3);
    entry = entry + Nx + 1;

    % B0(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j+1);
    vals(entries) = right_interface_A0(i,j,4);
    entry = entry + Nx + 1;

    for k = 1:Nx
        % An(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j-1) + k;
        vals(entries) = right_interface_An(k,:,i,j,1);
        entry = entry + Nx + 1;

        % Bn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j) + k;
        vals(entries) = right_interface_An(k,:,i,j,2);
        entry = entry + Nx + 1;

        % An(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j) + k;
        vals(entries) = right_interface_An(k,:,i,j,3);
        entry = entry + Nx + 1;

        % Bn(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j+1) + k;
        vals(entries) = right_interface_An(k,:,i,j,4);
        entry = entry + Nx + 1;
    end

    for k = 1:Ny
        % Cn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i-1,j) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,1);
        entry = entry + Nx + 1;

        % Cn(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i-1,j+1) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,3);
        entry = entry + Nx + 1;
    end

end


for i = 2:W-1

    % Left side blocks (not including corner blocks)
    j = 1;

    row_entries = row_index(i,j) + inc;

    % B0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j);
    vals(entries) = right_interface_A0(i,j,2);
    entry = entry + Nx + 1;

    % A0(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j);
    vals(entries) = right_interface_A0(i,j,3);
    entry = entry + Nx + 1;

    % B0(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j+1);
    vals(entries) = right_interface_A0(i,j,4);
    entry = entry + Nx + 1;

    for k = 1:Nx

        % Bn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j) + k;
        vals(entries) = right_interface_An(k,:,i,j,2);
        entry = entry + Nx + 1;

        % An(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j) + k;
        vals(entries) = right_interface_An(k,:,i,j,3);
        entry = entry + Nx + 1;

        % Bn(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j+1) + k;
        vals(entries) = right_interface_An(k,:,i,j,4);
        entry = entry + Nx + 1;

    end

    for k = 1:Ny

        % Cn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i-1,j) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,1);
        entry = entry + Nx + 1;

        % Dn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,2);
        entry = entry + Nx + 1;

        % Cn(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i-1,j+1) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,3);
        entry = entry + Nx + 1;

        % Dn(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j+1) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,4);
        entry = entry + Nx + 1;

    end


    % Right side blocks (not including corner blocks)
    j = L-1;

    row_entries = row_index(i,j) + inc;

    % A0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j-1);
    vals(entries) = right_interface_A0(i,j,1);
    entry = entry + Nx + 1;

    % B0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j);
    vals(entries) = right_interface_A0(i,j,2);
    entry = entry + Nx + 1;

    % A0(i,j+1) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j);
    vals(entries) = right_interface_A0(i,j,3);
    entry = entry + Nx + 1;

    for k = 1:Nx

        % An(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j-1) + k;
        vals(entries) = right_interface_An(k,:,i,j,1);
        entry = entry + Nx + 1;

        % Bn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j) + k;
        vals(entries) = right_interface_An(k,:,i,j,2);
        entry = entry + Nx + 1;

        % An(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j) + k;
        vals(entries) = right_interface_An(k,:,i,j,3);
        entry = entry + Nx + 1;

    end

    for k = 1:Ny

        % Cn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i-1,j) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,1);
        entry = entry + Nx + 1;

        % Dn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,2);
        entry = entry + Nx + 1;

        % Cn(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i-1,j+1) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,3);
        entry = entry + Nx + 1;

        % Dn(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j+1) + k;
        vals(entries) = right_interface_Bn(k,:,i,j,4);
        entry = entry + Nx + 1;

    end

    % Interior blocks (not including any boundary blocks)
    for j = 2:L-2

        row_entries = row_index(i,j) + inc;

        % A0(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j-1);
        vals(entries) = right_interface_A0(i,j,1);
        entry = entry + Nx + 1;

        % B0(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j);
        vals(entries) = right_interface_A0(i,j,2);
        entry = entry + Nx + 1;

        % A0(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j);
        vals(entries) = right_interface_A0(i,j,3);
        entry = entry + Nx + 1;

        % B0(i,j+1) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j+1);
        vals(entries) = right_interface_A0(i,j,4);
        entry = entry + Nx + 1;

        for k = 1:Nx

            % An(i,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_right(i,j-1) + k;
            vals(entries) = right_interface_An(k,:,i,j,1);
            entry = entry + Nx + 1;

            % Bn(i,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_right(i,j) + k;
            vals(entries) = right_interface_An(k,:,i,j,2);
            entry = entry + Nx + 1;

            % An(i,j+1) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_right(i,j) + k;
            vals(entries) = right_interface_An(k,:,i,j,3);
            entry = entry + Nx + 1;

            % Bn(i,j+1) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_right(i,j+1) + k;
            vals(entries) = right_interface_An(k,:,i,j,4);
            entry = entry + Nx + 1;

        end

        for k = 1:Ny

            % Cn(i,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_top(i-1,j) + k;
            vals(entries) = right_interface_Bn(k,:,i,j,1);
            entry = entry + Nx + 1;

            % Dn(i,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_top(i,j) + k;
            vals(entries) = right_interface_Bn(k,:,i,j,2);
            entry = entry + Nx + 1;

            % Cn(i,j+1) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_top(i-1,j+1) + k;
            vals(entries) = right_interface_Bn(k,:,i,j,3);
            entry = entry + Nx + 1;

            % Dn(i,j+1) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_top(i,j+1) + k;
            vals(entries) = right_interface_Bn(k,:,i,j,4);
            entry = entry + Nx + 1;

        end

    end

end



% Top-side boundary continuity conditions
row_index = @(i,j) (i - 1 + (j-1)*(W-1))*Ny + (L-1)*W*(Nx+1);
inc = 1:Ny;

% Bottom-left corner block
i = 1;
j = 1;

row_entries = row_index(i,j) + inc;

% B0(i,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j);
vals(entries) = top_interface_A0(:,i,j,2);
entry = entry + Ny;

% B0(i+1,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i+1,j);
vals(entries) = top_interface_A0(:,i,j,4);
entry = entry + Ny;

for k = 1:Nx

    % Bn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j) + k;
    vals(entries) = top_interface_An(k,:,i,j,2);
    entry = entry + Ny;

    % Bn(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i+1,j) + k;
    vals(entries) = top_interface_An(k,:,i,j,4);
    entry = entry + Ny;
end

for k = 1:Ny
    % Dn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,2);
    entry = entry + Ny;

    % Cn(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,3);
    entry = entry + Ny;

    % Dn(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i+1,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,4);
    entry = entry + Ny;
end


% Bottom-right corner block
i = 1;
j = L;

row_entries = row_index(i,j) + inc;

% A0(i,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j-1);
vals(entries) = top_interface_A0(:,i,j,1);
entry = entry + Ny;

% A0(i+1,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i+1,j-1);
vals(entries) = top_interface_A0(:,i,j,3);
entry = entry + Ny;

for k = 1:Nx
    % An(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j-1) + k;
    vals(entries) = top_interface_An(k,:,i,j,1);
    entry = entry + Ny;

    % An(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i+1,j-1) + k;
    vals(entries) = top_interface_An(k,:,i,j,3);
    entry = entry + Ny;
end

for k = 1:Ny
    % Dn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,2);
    entry = entry + Ny;

    % Cn(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,3);
    entry = entry + Ny;

    % Dn(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i+1,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,4);
    entry = entry + Ny;
end


% Top-left corner block
i = W-1;
j = 1;

row_entries = row_index(i,j) + inc;

% B0(i,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j);
vals(entries) = top_interface_A0(:,i,j,2);
entry = entry + Ny;

% B0(i+1,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i+1,j);
vals(entries) = top_interface_A0(:,i,j,4);
entry = entry + Ny;

for k = 1:Nx
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j) + k;
    vals(entries) = top_interface_An(k,:,i,j,2);
    entry = entry + Ny;

    % Bn(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i+1,j) + k;
    vals(entries) = top_interface_An(k,:,i,j,4);
    entry = entry + Ny;
end

for k = 1:Ny
    % Cn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i-1,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,1);
    entry = entry + Ny;

    % Dn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,2);
    entry = entry + Ny;

    % Cn(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,3);
    entry = entry + Ny;
end


% Top-right corner block
i = W-1;
j = L;

row_entries = row_index(i,j) + inc;

% A0(i,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i,j-1);
vals(entries) = top_interface_A0(:,i,j,1);
entry = entry + Ny;

% A0(i+1,j) values
entries = inc + entry;
rows(entries) = row_entries;
cols(entries) = index_right(i+1,j-1);
vals(entries) = top_interface_A0(:,i,j,3);
entry = entry + Ny;

for k = 1:Nx
    % An(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j-1) + k;
    vals(entries) = top_interface_An(k,:,i,j,1);
    entry = entry + Ny;

    % An(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i+1,j-1) + k;
    vals(entries) = top_interface_An(k,:,i,j,3);
    entry = entry + Ny;
end

for k = 1:Ny
    % Cn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i-1,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,1);
    entry = entry + Ny;

    % Dn(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,2);
    entry = entry + Ny;

    % Cn(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_top(i,j) + k;
    vals(entries) = top_interface_Bn(k,:,i,j,3);
    entry = entry + Ny;
end


for i = 2:W-2
    % Left side blocks (not including corner blocks)
    j = 1;

    row_entries = row_index(i,j) + inc;

    % B0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j);
    vals(entries) = top_interface_A0(:,i,j,2);
    entry = entry + Ny;

    % B0(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i+1,j);
    vals(entries) = top_interface_A0(:,i,j,4);
    entry = entry + Ny;

    for k = 1:Nx
        % Bn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j) + k;
        vals(entries) = top_interface_An(k,:,i,j,2);
        entry = entry + Ny;

        % Bn(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i+1,j) + k;
        vals(entries) = top_interface_An(k,:,i,j,4);
        entry = entry + Ny;
    end

    for k = 1:Ny
        % Cn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i-1,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,1);
        entry = entry + Ny;

        % Dn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,2);
        entry = entry + Ny;

        % Cn(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,3);
        entry = entry + Ny;

        % Dn(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i+1,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,4);
        entry = entry + Ny;
    end


    % Right side blocks (not including corner blocks)
    j = L;

    row_entries = row_index(i,j) + inc;

    % A0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j-1);
    vals(entries) = top_interface_A0(:,i,j,1);
    entry = entry + Ny;

    % A0(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i+1,j-1);
    vals(entries) = top_interface_A0(:,i,j,3);
    entry = entry + Ny;

    for k = 1:Nx
        % An(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j-1) + k;
        vals(entries) = top_interface_An(k,:,i,j,1);
        entry = entry + Ny;

        % An(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i+1,j-1) + k;
        vals(entries) = top_interface_An(k,:,i,j,3);
        entry = entry + Ny;
    end

    for k = 1:Ny
        % Cn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i-1,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,1);
        entry = entry + Ny;

        % Dn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,2);
        entry = entry + Ny;

        % Cn(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,3);
        entry = entry + Ny;

        % Dn(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i+1,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,4);
        entry = entry + Ny;
    end
end


for j = 2:L-1
    % Bottom side blocks (not including corner blocks)
    i = 1;

    row_entries = row_index(i,j) + inc;

    % A0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j-1);
    vals(entries) = top_interface_A0(:,i,j,1);
    entry = entry + Ny;

    % B0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j);
    vals(entries) = top_interface_A0(:,i,j,2);
    entry = entry + Ny;

    % A0(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i+1,j-1);
    vals(entries) = top_interface_A0(:,i,j,3);
    entry = entry + Ny;

    % B0(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i+1,j);
    vals(entries) = top_interface_A0(:,i,j,4);
    entry = entry + Ny;

    for k = 1:Nx
        % An(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j-1) + k;
        vals(entries) = top_interface_An(k,:,i,j,1);
        entry = entry + Ny;

        % Bn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j) + k;
        vals(entries) = top_interface_An(k,:,i,j,2);
        entry = entry + Ny;

        % An(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i+1,j-1) + k;
        vals(entries) = top_interface_An(k,:,i,j,3);
        entry = entry + Ny;

        % Bn(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i+1,j) + k;
        vals(entries) = top_interface_An(k,:,i,j,4);
        entry = entry + Ny;
    end

    for k = 1:Ny
        % Dn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,2);
        entry = entry + Ny;

        % Cn(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,3);
        entry = entry + Ny;

        % Dn(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i+1,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,4);
        entry = entry + Ny;
    end


    % Top side blocks (not including corner blocks)
    i = W-1;

    row_entries = row_index(i,j) + inc;

    % A0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j-1);
    vals(entries) = top_interface_A0(:,i,j,1);
    entry = entry + Ny;

    % B0(i,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i,j);
    vals(entries) = top_interface_A0(:,i,j,2);
    entry = entry + Ny;

    % A0(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i+1,j-1);
    vals(entries) = top_interface_A0(:,i,j,3);
    entry = entry + Ny;

    % B0(i+1,j) values
    entries = inc + entry;
    rows(entries) = row_entries;
    cols(entries) = index_right(i+1,j);
    vals(entries) = top_interface_A0(:,i,j,4);
    entry = entry + Ny;

    for k = 1:Nx
        % An(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j-1) + k;
        vals(entries) = top_interface_An(k,:,i,j,1);
        entry = entry + Ny;

        % Bn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j) + k;
        vals(entries) = top_interface_An(k,:,i,j,2);
        entry = entry + Ny;

        % An(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i+1,j-1) + k;
        vals(entries) = top_interface_An(k,:,i,j,3);
        entry = entry + Ny;

        % Bn(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i+1,j) + k;
        vals(entries) = top_interface_An(k,:,i,j,4);
        entry = entry + Ny;
    end

    for k = 1:Ny
        % Cn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i-1,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,1);
        entry = entry + Ny;

        % Dn(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,2);
        entry = entry + Ny;

        % Cn(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_top(i,j) + k;
        vals(entries) = top_interface_Bn(k,:,i,j,3);
        entry = entry + Ny;
    end


    % Interior blocks (not including any boundary blocks)
    for i = 2:W-2
        row_entries = row_index(i,j) + inc;

        % A0(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j-1);
        vals(entries) = top_interface_A0(:,i,j,1);
        entry = entry + Ny;

        % B0(i,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i,j);
        vals(entries) = top_interface_A0(:,i,j,2);
        entry = entry + Ny;

        % A0(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i+1,j-1);
        vals(entries) = top_interface_A0(:,i,j,3);
        entry = entry + Ny;

        % B0(i+1,j) values
        entries = inc + entry;
        rows(entries) = row_entries;
        cols(entries) = index_right(i+1,j);
        vals(entries) = top_interface_A0(:,i,j,4);
        entry = entry + Ny;

        for k = 1:Nx
            % An(i,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_right(i,j-1) + k;
            vals(entries) = top_interface_An(k,:,i,j,1);
            entry = entry + Ny;

            % Bn(i,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_right(i,j) + k;
            vals(entries) = top_interface_An(k,:,i,j,2);
            entry = entry + Ny;

            % An(i+1,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_right(i+1,j-1) + k;
            vals(entries) = top_interface_An(k,:,i,j,3);
            entry = entry + Ny;

            % Bn(i+1,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_right(i+1,j) + k;
            vals(entries) = top_interface_An(k,:,i,j,4);
            entry = entry + Ny;
        end

        for k = 1:Ny
            % Cn(i,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_top(i-1,j) + k;
            vals(entries) = top_interface_Bn(k,:,i,j,1);
            entry = entry + Ny;

            % Dn(i,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_top(i,j) + k;
            vals(entries) = top_interface_Bn(k,:,i,j,2);
            entry = entry + Ny;

            % Cn(i+1,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_top(i,j) + k;
            vals(entries) = top_interface_Bn(k,:,i,j,3);
            entry = entry + Ny;

            % Dn(i+1,j) values
            entries = inc + entry;
            rows(entries) = row_entries;
            cols(entries) = index_top(i+1,j) + k;
            vals(entries) = top_interface_Bn(k,:,i,j,4);
            entry = entry + Ny;
        end
    end
end

A = sparse(rows,cols,vals,grid.total_unknowns,grid.total_unknowns,total_entries);

end