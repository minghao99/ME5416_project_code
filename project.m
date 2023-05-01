home_theta_list=[0.2;0.4;0.5;0.6;0.8;1.2;1.5];

body1 = rigidBody('body1');
jnt1 = rigidBodyJoint('jnt1','revolute');
jnt1.HomePosition = home_theta_list(1);
jnt1.JointAxis = [0,0,1];
tform = trvec2tform([0, 0, 0.17]);
setFixedTransform(jnt1,tform);
body1.Joint = jnt1;
robot = rigidBodyTree;
addBody(robot,body1,'base')

body2 = rigidBody('body2');
jnt2 = rigidBodyJoint('jnt2','revolute');
jnt2.HomePosition = home_theta_list(2);
jnt2.JointAxis = [1,0,0];
tform2 = trvec2tform([0, 0, 0.17]);
setFixedTransform(jnt2,tform2);
body2.Joint = jnt2;
addBody(robot,body2,'body1');

body3 = rigidBody('body3');
jnt3 = rigidBodyJoint('jnt3','revolute');
tform3 = trvec2tform([0, 0, 0.2]);
setFixedTransform(jnt3,tform3);
jnt3.HomePosition = home_theta_list(3);
body3.Joint = jnt3;
jnt3.JointAxis = [0,0,1];
addBody(robot,body3,'body2');

body4 = rigidBody('body4');
jnt4 = rigidBodyJoint('jnt4','revolute');
tform4 = trvec2tform([0, 0, 0.2]);
setFixedTransform(jnt4,tform4);
jnt4.HomePosition = home_theta_list(4);
jnt4.JointAxis = [1,0,0];
body4.Joint = jnt4;
addBody(robot,body4,'body3');

body5 = rigidBody('body5');
jnt5 = rigidBodyJoint('jnt5','revolute');
tform5 = trvec2tform([0, 0, 0.2]);
setFixedTransform(jnt5,tform5);
jnt5.HomePosition = home_theta_list(5);
jnt5.JointAxis = [0,0,1];
body5.Joint = jnt5;
addBody(robot,body5,'body4');

body6 = rigidBody('body6');
jnt6 = rigidBodyJoint('jnt6','revolute');
tform6 = trvec2tform([0, 0, 0.2]);
setFixedTransform(jnt6,tform6);
jnt6.HomePosition = home_theta_list(6);
jnt6.JointAxis = [1,0,0];
body6.Joint = jnt6;
addBody(robot,body6,'body5');

body7 = rigidBody('body7');
jnt7 = rigidBodyJoint('jnt7','revolute');
tform7 = trvec2tform([0, 0, 0.126]);
setFixedTransform(jnt7,tform7);
jnt7.HomePosition = home_theta_list(7);
jnt7.JointAxis = [0,0,1];
body7.Joint = jnt7;
addBody(robot,body7,'body6');

bodyEndEffector = rigidBody('endeffector');
tform8 = trvec2tform([0, 0, 0.1]);
setFixedTransform(bodyEndEffector.Joint,tform8);
addBody(robot,bodyEndEffector,'body7');

viztree = interactiveRigidBodyTree(robot);

all_screw_axis = transpose([0 0 1 0 0 0;
                  1 0 0 0 0.34 0;
                  0 0 1 0 0 0;
                  1 0 0 0 0.74 0;
                  0 0 1 0 0 0;
                  1 0 0 0 1.14 0;
                  0 0 1 0 0 0]);
% theta_list = [0;0;0;0;0;0;0];
zero_T = [eye(3) [0;0;1.366]; 0 0 0 1];

Theta_0 = [homeConfiguration(robot).JointPosition]';
space_J0 = space_jacobian(all_screw_axis,Theta_0);
space_T0 = space_forward_kinematics(all_screw_axis,Theta_0,zero_T);
body_J0 = body_jacobian(space_T0,space_J0);

Theta_cur = Theta_0;
space_T_cur = space_T0;
space_J_cur = space_J0;
body_J_cur = body_J0;

theta_test = [0.2;0.4;0.5;0.6;0.8;1.2;1.5] + [0.5;-0.5;0.5;0.5;0.5;-0.5;0.5];
space_T_test = space_forward_kinematics(all_screw_axis,theta_test,zero_T);
space_Goal_T = space_T_test;

% space_Goal_T = [[cos(pi/3) -sin(pi/3) 0; sin(pi/3) cos(pi/3) 0; 0 0 1] [0;0;1.366]; 0 0 0 1]; % user define

% T0\space_Goal_T is the same as inverse(T0)*space_Goal_T

% space V = T_dot * T_inverse
space_V0_matrix = (space_Goal_T - space_T_cur)*inv(space_T_cur);
space_V_cur = twist_matrix_to_R6(space_V0_matrix);



step = 0.05;
for t = 0:step:10
    % Compute next step joint config using current inverse jacobian mat
    % based on current joint config. left Pseudo inverse, as
    % underconstrained. pinv(J) = (J'*J)^-1 * J'

    Theta_cur = Theta_cur + step*pinv(space_J_cur,0.05)*space_V_cur;

    % Update jacobian mat for next step
    space_J_cur = space_jacobian(all_screw_axis,Theta_cur);


    % Move the displayed robot arm to next joint config
    viztree.Configuration = Theta_cur;
    
    % Update EE position and re-calculate intended EE velocity
    space_T_cur = space_forward_kinematics(all_screw_axis,Theta_cur,zero_T);
    space_V_cur_matrix = (space_Goal_T - space_T_cur)*inv(space_T_cur);
    space_V_cur = twist_matrix_to_R6(space_V_cur_matrix);
    
    speed = size_twist(space_V_cur);

    space_V_cur=space_V_cur/speed;

    if speed < 0.01
        break
    end
    % Pause to simulate closer to real time
    pause(0.1)
end

















function R3_matrix = R3_to_matrix(W)
    R3_matrix = [0 -W(3) W(2); W(3) 0 -W(1); -W(2) W(1) 0];
end

function exp_of_screw = exp_screw(screw_axis,theta)
    omega = screw_axis(1:3);
    v = screw_axis(4:6);
    if sqrt(sum(omega.*omega)) == 1
        omega_mat = R3_to_matrix(screw_axis(1:3));
        G_theta = eye(3)*theta + (1-cos(theta))*omega_mat + (theta - sin(theta))*omega_mat*omega_mat;
        exp_of_screw = [expm(omega_mat*theta) G_theta*v; 0 0 0 1];
    end

    if sqrt(sum(omega.*omega)) == 0 && sqrt(sum(v.*v)) == 1
        exp_of_screw = [eye(3) theta*v; 0 0 0 1];
    end
end

% adjoint_T is the 6x6 adjoint matrix of 4x4 homogeneous transformation
%   matrix T
function adjoint_T = adjoint(T)
    R = T(1:3,1:3);
    P = T(1:3,4);
    adjoint_T = [R eye(3)*0;
                 R3_to_matrix(P)*R R];
end

% S_list is 6xn, where each column is the screw axis for each link
% theta_list = nx1, each row is the joint displacement from 'zero' position
% zero_T is the EE transformation when at 'zero' position
function space_forward_T = space_forward_kinematics(S_list,theta_list,zero_T)
    n = length(theta_list);
    space_forward_T = eye(4);
    for k = 1:n
        space_forward_T = space_forward_T * exp_screw(S_list(:,k),theta_list(k));
    end
    space_forward_T = space_forward_T * zero_T;
end

% space_J is 6xn matrix, multiplying with nx1 theta_dot_vector to give
%  space 6x1 twist of EE.
function space_J = space_jacobian(S_list,theta_list)
    n = length(theta_list);
    J1 = S_list(:,1);
    space_J = J1;
    for k = 2:n
        T = eye(4);
        for i = 1:(k-1)
            T = T*exp_screw(S_list(:,i),theta_list(i));
        end
        space_J = [space_J adjoint(T)*S_list(:,k)];
    end
end

% Body jacobian can be easily obtained from space jacobian
% body jacobian = adjoint(inverse of space_T_end_effector) * space jacobian
function body_J = body_jacobian(forward_T,space_J)
    body_J = adjoint(inv(forward_T))*space_J;
end

function R6_twist = twist_matrix_to_R6(T)
    temp = T;
    R6_twist = [temp(3,2);temp(1,3);temp(2,1);temp(1,4);temp(2,4);temp(3,4)];
end

function size_of_twist = size_twist(twist)
    omega = twist(1:3);
    v = twist(4:6);
    if sqrt(sum(omega.*omega)) == 0
        size_of_twist = sqrt(sum(v.*v));
    else
        size_of_twist = sqrt(sum(omega.*omega));
    end
end