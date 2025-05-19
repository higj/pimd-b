#pragma once

class WalkersCommunicationBase {
public:
    explicit WalkersCommunicationBase();
    virtual ~WalkersCommunicationBase() = default;
    
    virtual void communicate();
};