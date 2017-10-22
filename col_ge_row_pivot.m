function [order,A] = col_ge_row_pivot(A)

    [~,n] = size(A);
    
    order = zeros(n,1);
    original_order = (1:n);
    
    % For each column
    for i = 1:n
        % Find pivot in row
        [~, index] = max(abs(A(:,i)));
        order(i) = original_order(index);
        % Swap top-most row with pivot row
        if (index~=i)
            temp = A(i,:);
            A(i,:) = A(index,:);
            A(index,:) = temp;
            
            temp = original_order(i);
            original_order(i) = original_order(index);
            original_order(index) = temp;
        end
        % Perform elimination
        for j = i+1:n
            A(:,j) = A(:,j) - (A(i,j)/A(i,i))*A(:,i);
        end
    end

end