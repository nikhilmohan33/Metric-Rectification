function rectifyImage(filename,debug)
    img = imread(filename);
    [d,fname] = fileparts(filename);
    rectI = rectifyAffine(filename,debug);
    affine_out = strcat('images/',fname,'_affine','.jpg');
    imshow(rectI);
    imwrite(rectI,affine_out);
    final_img = rectifyMetric(affine_out,debug);
    figure;
    imshow(final_img)
    metric_out = strcat('images/',fname,'_metric','.jpg');
    imwrite(final_img,metric_out);
end

function rectI = rectifyAffine(filename,debug)
    img = imread(filename);
    [d,fname] = fileparts(filename);
    if(debug)
        figure;
        imshow(img);
        hold on;
        fullfile('images',[fname '.mat'])
        [x y] = ginput(8); 
        save(strcat('data/',fname),'x','y');
        for i = 1:4
            if(i > 2)
                color = 'r'
            else
                color = 'b'
            end
            line([x(2*i-1),x(2*i)],[y(2*i-1),y(2*i)],'color',color,'LineWidth',2)
        end
        saveas(gcf,strcat('images/',fname,'_parallel','.jpg'))
    end
    data = load(strcat('data/',fname,'.mat'));
    x = data.x;
    y = data.y;
    l1 = cross([x(1) y(1) 1],[x(2) y(2) 1]);
    l2 = cross([x(3) y(3) 1],[x(4) y(4) 1]);
    l3 = cross([x(5) y(5) 1],[x(6) y(6) 1]);
    l4 = cross([x(7) y(7) 1],[x(8) y(8) 1]);
    p1 = cross(l1,l2);
    p1 = p1/p1(3);
    p2 = cross(l3,l4);
    p2 = p2/p2(3);
    l_inf = cross(p1,p2);
    figure;
    imshow(img);
    line([p1(1) p2(1)], [p1(2) p2(2)],'LineWidth',2);
    axis auto
    saveas(gcf,strcat('images/',fname,'_vanishing_line','.jpg'))
    H = eye(3);
    H(3,1) = l_inf(1)/l_inf(3);
    H(3,2) = l_inf(2)/l_inf(3);
    H(3,3) = 1;
    rectI = applyH(img,H);
end


function final_img = rectifyMetric(filename,debug)
    img = imread(filename);
    [d,fname] = fileparts(filename);
    if(debug)
        figure;
        imshow(img);
        hold on;
        fullfile('images',[fname '.mat'])
        [x y] = ginput(8); 
        save(strcat('data/',fname),'x','y');
        for i = 1:4
            if(i < 3)
                color = 'r'
            else
                color = 'g'
            end
            line([x(2*i-1),x(2*i)],[y(2*i-1),y(2*i)],'color',color,'LineWidth',2)
        end
        saveas(gcf,strcat('images/',fname,'_orthogonal','.jpg'))
    end
    data = load(strcat('data/',fname,'.mat'));
    x = data.x;
    y = data.y;
    l1 = cross([x(1) y(1) 1],[x(2) y(2) 1]);
    m1 = cross([x(3) y(3) 1],[x(4) y(4) 1]);
    l2 = cross([x(5) y(5) 1],[x(6) y(6) 1]);
    m2 = cross([x(7) y(7) 1],[x(8) y(8) 1]);
    
    cosTheta1 = dot([l1(1) l1(2)],[m1(1) m1(2)])/(norm([l1(1) l1(2)])*norm([m1(1) m1(2)]));
    cosTheta2 = dot([l2(1) l2(2)],[m2(1) m2(2)])/(norm([l2(1) l2(2)])*norm([m2(1) m2(2)]));    
    before = [cosTheta1 cosTheta2];
    l11 = [l1(1)/l1(3) l1(2)/l1(3)];
    m11 = [m1(1)/m1(3) m1(2)/m1(3)];
    l22 = [l2(1)/l2(3) l2(2)/l2(3)];
    m22 = [m2(1)/m2(3) m2(2)/m2(3)];
    
    A = [l11(1)*m11(1),l11(1)*m11(2)+l11(2)*m11(1);
        l22(1)*m22(1),l22(1)*m22(2)+l22(2)*m22(1)];
    
    b = [-l11(2)*m11(2);-l22(2)*m22(2)];
    
    x = A \ b;
    
    S = [x(1) x(2); x(2) 1];
    [U,D,V] = svd(S);
    
    A = U*sqrt(D)*V';
    H = eye(3);
    H(1,1) = A(1,1);
    H(1,2) = A(1,2);
    H(2,1) = A(2,1);
    H(2,2) = A(2,2);
    
    l1_hat = H * l1';
    m1_hat = H * m1';
    l2_hat = H * l2';
    m2_hat = H * m2';
    newcos1 = dot([l1_hat(1) l1_hat(2)],[m1_hat(1) m1_hat(2)])/(norm([l1(1) l1(2)])*norm([m1(1) m1(2)]));
    newcos2 = dot([l2_hat(1) l2_hat(2)],[m2_hat(1) m2_hat(2)])/(norm([l2(1) l2(2)])*norm([m2(1) m2(2)]));
    after = [newcos1 newcos2];
    
    
    
    before
    after
    final_img = applyH(img,inv(H));
end
