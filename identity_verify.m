function accept = identity_verify(identity_tx, identity_rx, q, k, delta, j)
%IDENTITY_VERIFY Verify one identity against one transmitted ID tag.
%
%   accept = identity_verify(identity_tx,identity_rx,q,k,delta,j) computes
%   the transmitted tag T_identity_tx(j) and the receiver-side verification
%   tag T_identity_rx(j). The receiver accepts exactly when the two tags are
%   equal.

    transmitted_tag = concatenated_rs_tag(identity_tx, q, k, delta, j);
    verification_tag = concatenated_rs_tag(identity_rx, q, k, delta, j);
    accept = (transmitted_tag == verification_tag);
end
