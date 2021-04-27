
function solution_error(num_pts, u_extended)
    x = range(0,1,length=num_pts)
    y = range(0,1,length=num_pts)'

    u_true = zeros((num_pts, num_pts));

    for ix = 1:num_pts
        for iy = 1:num_pts
            u_true[ix,iy]=u(x[ix],y[iy])
        end
    end

    error = u_true-u_extended;
    error = error[:]
    return x,y, error
end

