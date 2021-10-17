function [p, g, l, l_dot, r, r_dot, aa_new, ad_new, am_new, S, phi, phi_dot] = RK4_new(p, g, l, l_dot, r, r_dot, a_aa, h, c, M, S, V)
    aa = -3*l_dot(1)*r_dot(1);
    am = -5*l_dot(1)*r_dot(1);
    ad = -(2*r_dot(2) - c*r(2))*l_dot(2)/cos(g(2) - l(2))  +  (r(2)/r(1))*(2*r_dot(1) - c*r(1))*l_dot(1)/cos(g(2)-l(2))  +  (r(2)/r(1))*cos(g(1) - l(1))*aa/cos(g(2)-l(2))  -  M*sign(S)/cos(g(2)-l(2)); 
    
    a_a = [aa ; ad; am];
    a = [aa ;  ad ;  am];
    
    g_b = g + (h/2)*[a(1)/V(1); a(2)/V(2) ; a(3)/V(3)];
    p_b = p + (h/2)*[V(1)*cos(g_b(1)) , V(1)*sin(g_b(1)) ; V(2)*cos(g_b(2)), V(2)*sin(g_b(2)) ; V(3)*cos(g_b(3)) , V(3)*sin(g_b(3)) ];
    l = [atan((p_b(3,2)-p_b(1,2))/(p_b(3,1)-p_b(1,1))) ;  atan((p_b(3,2)-p_b(2,2))/(p_b(3,1)-p_b(2,1)))];
    r = [((p_b(1,1)-p_b(3,1))^2 + (p_b(1,2)-p_b(3,2))^2)^0.5 ; ((p_b(2,1)-p_b(3,1))^2 + (p_b(2,2)-p_b(3,2))^2)^0.5];
    r_dot = [V(3)*cos(g_b(3)-l(1)) - V(1)*cos(g_b(1)-l(1)) ; V(3)*cos(g_b(3) - l(2)) - V(2)*cos(g_b(2) - l(2))];
    l_dot = [ ( V(3)*sin(g_b(3) - l(1)) - V(1)*sin(g_b(1) - l(1)) )/r(1) ; ( V(3)*sin(g_b(3) - l(2)) - V(2)*sin(g_b(2) - l(2)) )/r(2)];
    %S = l_dot(1) - l_dot(2) + c*(l(1) - l(2));
    
    a(1) = -3*l_dot(1)*r_dot(1);
    a(3) = -5*l_dot(1)*r_dot(1);
    a(2) = -(2*r_dot(2) - c*r(2))*l_dot(2)/cos(g_b(2) - l(2))  +  (r(2)/r(1))*(2*r_dot(1) - c*r(1))*l_dot(1)/cos(g_b(2)-l(2))  +  (r(2)/r(1))*cos(g_b(1) - l(1))*a(1)/cos(g_b(2)-l(2))  -  M*sign(S)/cos(g_b(2)-l(2)) ; 
    
    a_b = [a(1); a(2) ; a(3)];
    
    g_c = g + (h/2)*[a(1)/V(1); a(2)/V(2) ; a(3)/V(3)];
    p_c = p + (h/2)*[V(1)*cos(g_c(1)) , V(1)*sin(g_c(1)) ; V(2)*cos(g_c(2)) ,  V(2)*sin(g_c(2)) ; V(3)*cos(g_c(3)) , V(3)*sin(g_c(3)) ];
    l = [atan((p_c(3,2)-p_c(1,2))/(p_c(3,1)-p_c(1,1))) ;  atan((p_c(3,2)-p_c(2,2))/(p_c(3,1)-p_c(2,1)))];
    r = [((p_c(1,1)-p_c(3,1))^2 + (p_c(1,2)-p_c(3,2))^2)^0.5 ; ((p_c(2,1)-p_c(3,1))^2 + (p_c(2,2)-p_c(3,2))^2)^0.5];
    r_dot = [V(3)*cos(g_c(3)-l(1)) - V(1)*cos(g_c(1)-l(1)) ; V(3)*cos(g_c(3) - l(2)) - V(2)*cos(g_c(2) - l(2))];
    l_dot = [ ( V(3)*sin(g_c(3) - l(1)) - V(1)*sin(g_c(1) - l(1)) )/r(1) ; ( V(3)*sin(g_c(3) - l(2)) - V(2)*sin(g_c(2) - l(2)) )/r(2)];
    %S = l_dot(1) - l_dot(2) + c*(l(1) - l(2));
    
    a(1) = -3*l_dot(1)*r_dot(1);
    a(3) = -5*l_dot(1)*r_dot(1);
    a(2) = -(2*r_dot(2) - c*r(2))*l_dot(2)/cos(g_c(2) - l(2))  +  (r(2)/r(1))*(2*r_dot(1) - c*r(1))*l_dot(1)/cos(g_c(2)-l(2))  +  (r(2)/r(1))*cos(g_c(1) - l(1))*a(1)/cos(g_c(2)-l(2))  -  M*sign(S)/cos(g_c(2)-l(2)); 
    
    a_c = [a(1); a(2) ; a(3)];
    
    
    g_d = g + h*[a(1)/V(1); a(2)/V(2) ; a(3)/V(3)];
    p_d = p + h*[V(1)*cos(g_d(1)) , V(1)*sin(g_d(1)) ; V(2)*cos(g_d(2)) ,  V(2)*sin(g_d(2)) ; V(3)*cos(g_d(3)) , V(3)*sin(g_d(3)) ];
    l = [atan((p_d(3,2)-p_d(1,2))/(p_d(3,1)-p_d(1,1))) ;  atan((p_d(3,2)-p_d(2,2))/(p_d(3,1)-p_d(2,1)))];
    r = [((p_d(1,1)-p_d(3,1))^2 + (p_d(1,2)-p_d(3,2))^2)^0.5 ; ((p_d(2,1)-p_d(3,1))^2 + (p_d(2,2)-p_d(3,2))^2)^0.5];
    r_dot = [V(3)*cos(g_d(3)-l(1)) - V(1)*cos(g_d(1)-l(1)) ; V(3)*cos(g_d(3) - l(2)) - V(2)*cos(g_d(2) - l(2))];
    l_dot = [ ( V(3)*sin(g_d(3) - l(1)) - V(1)*sin(g_d(1) - l(1)) )/r(1) ; ( V(3)*sin(g_d(3) - l(2)) - V(2)*sin(g_d(2) - l(2)) )/r(2)];
    %S = l_dot(1) - l_dot(2) + c*(l(1) - l(2));
    
    a(1) = -3*l_dot(1)*r_dot(1);
    a(3) = -5*l_dot(1)*r_dot(1);
    a(2) = -(2*r_dot(2) - c*r(2))*l_dot(2)/cos(g_d(2) - l(2))  +  (r(2)/r(1))*(2*r_dot(1) - c*r(1))*l_dot(1)/cos(g_d(2)-l(2))  +  (r(2)/r(1))*cos(g_d(1) - l(1))*a(1)/cos(g_d(2)-l(2))  -  M*sign(S)/cos(g_d(2)-l(2)) ; 
    
    a_d = [a(1); a(2) ; a(3)];
    
    %Final values to be returned
    go = g + (h/6)*( [a_a(1)/V(1) ; a_a(2)/V(2) ; a_a(3)/V(3)] + 2*[a_b(1)/V(1) ; a_b(2)/V(2) ; a_b(3)/V(3)] + 2*[a_c(1)/V(1) ; a_c(2)/V(2) ; a_c(3)/V(3)] + [a_d(1)/V(1) ; a_d(2)/V(2) ; a_d(3)/V(3)]);
    p = p + (h/6)*( [V(1)*cos(g(1)) , V(1)*sin(g(1)) ; V(2)*cos(g(2)), V(2)*sin(g(2)) ; V(3)*cos(g(3)) , V(3)*sin(g(3))]+ 2*[V(1)*cos(g_b(1)) , V(1)*sin(g_b(1)) ; V(2)*cos(g_b(2)), V(2)*sin(g_b(2)) ; V(3)*cos(g_b(3)) , V(3)*sin(g_b(3)) ] + 2*[V(1)*cos(g_c(1)) , V(1)*sin(g_c(1)) ; V(2)*cos(g_c(2)), V(2)*sin(g_c(2)) ; V(3)*cos(g_c(3)) , V(3)*sin(g_c(3))]  + [V(1)*cos(g_d(1)) , V(1)*sin(g_d(1)) ; V(2)*cos(g_d(2)), V(2)*sin(g_d(2)) ; V(3)*cos(g_d(3)) , V(3)*sin(g_d(3))] );
    g = go;
    
    l = [atan((p(3,2)-p(1,2))/(p(3,1)-p(1,1))) ;  atan((p(3,2)-p(2,2))/(p(3,1)-p(2,1)))];
    r = [((p(1,1)-p(3,1))^2 + (p(1,2)-p(3,2))^2)^0.5 ; ((p(2,1)-p(3,1))^2 + (p(2,2)-p(3,2))^2)^0.5];
    r_dot = [V(3)*cos(g(3)-l(1)) - V(1)*cos(g(1)-l(1)) ; V(3)*cos(g(3) - l(2)) - V(2)*cos(g(2) - l(2))];
    l_dot = [ ( V(3)*sin(g(3) - l(1)) - V(1)*sin(g(1) - l(1)) )/r(1) ; ( V(3)*sin(g(3) - l(2)) - V(2)*sin(g(2) - l(2)) )/r(2)];
    
    phi_dot = l_dot(1) - l_dot(2);
    phi = l(1) - l(2);
    S = phi_dot + c*phi;
    
    aa_new = -3*l_dot(1)*r_dot(1);
    am_new = -5*l_dot(1)*r_dot(1);
    ad_new = -(2*r_dot(2) - c*r(2))*l_dot(2)/cos(g_d(2) - l(2))  +  (r(2)/r(1))*(2*r_dot(1) - c*r(1))*l_dot(1)/cos(g_d(2)-l(2))  +  (r(2)/r(1))*cos(g_d(1) - l(1))*aa_new/cos(g_d(2)-l(2))  -  M*sign(S)/cos(g_d(2)-l(2)) ; 
   
end
